"""Exercise concurrent scheduling without submitting or changing Slurm jobs."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock

spec = importlib.util.spec_from_file_location("overlap", Path(__file__).with_name("manage-study-overlap.py"))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class OverlapTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="imr-overlap-test-")
        self.addCleanup(self.temporary.cleanup)
        self.study = Path(self.temporary.name)
        self.c = module.Controller.__new__(module.Controller)
        self.c.study = self.study
        self.c.directory = self.study
        self.c.policy = {"pilot_ids": [1, 2], "full_ids": [3, 4], "concurrency": 32, "pilot_reserve": 2}
        self.c.state = {"retried": {}, "submission_in_progress": None, "blocked": None}
        for phase, job, ids in [("pilot", "100", [1, 2]), ("full", "200", [3, 4])]:
            self.c.state[phase] = {"arrays": [dict(job=job, ids=ids, memory_mib=16384, hours=48)],
                                   "accepted": False, "validator": None, "throttle": 30}
        for name in ("save", "event", "retry", "validate", "verify", "route", "block", "command"):
            setattr(self.c, name, Mock())

    def complete(self, *ids):
        for task in ids:
            p = self.study / "tasks" / ("%04d" % task)
            p.mkdir(parents=True, exist_ok=True)
            (p / "TASK-COMPLETED").write_text("fixture\n")

    def test_concurrency_reserves_pilot_slots(self):
        self.assertEqual(self.c.production_slots(), 30)
        self.assertEqual(self.c.production_slots() + self.c.policy["pilot_reserve"], 32)
        self.c.state["pilot"]["accepted"] = True
        self.assertEqual(self.c.production_slots(), 32)

    def test_resource_retry_limits(self):
        self.assertEqual(module.retry_request({1: "OUT_OF_MEMORY"}, 16384, 48, {}), (32768, 48))
        self.assertEqual(module.retry_request({1: "TIMEOUT"}, 16384, 48, {}), (16384, 96))
        self.assertEqual(module.retry_request({1: "NODE_FAIL"}, 16384, 48, {}), (16384, 48))
        for failure, memory, hours, retried in [
            ("FAILED", 16384, 48, {}), ("CANCELLED", 16384, 48, {}),
            ("TIMEOUT", 16384, 96, {}), ("OUT_OF_MEMORY", 200000, 48, {}),
            ("NODE_FAIL", 16384, 48, {"1": {}})]:
            with self.subTest(failure=failure, hours=hours, retried=retried):
                with self.assertRaises(ValueError):
                    module.retry_request({1: failure}, memory, hours, retried)

    def test_running_pilot_does_not_block_production_computation(self):
        self.complete(2, 3, 4)
        self.c.advance("pilot", [["100_1", "RUNNING", "highmem", "48:00:00"]], {})
        self.c.advance("full", [], {"200_3": "COMPLETED", "200_4": "COMPLETED"})
        self.c.validate.assert_not_called()  # Final acceptance still waits for pilot.
        self.c.retry.assert_not_called()

    def test_pilot_validator_waits_for_terminal_array(self):
        self.complete(1, 2)
        self.c.advance("pilot", [["100_1", "COMPLETING", "highmem", "48:00:00"]], {})
        self.c.validate.assert_not_called()
        self.c.advance("pilot", [], {})
        self.c.validate.assert_called_once_with("pilot")

    def test_resource_failure_retries_only_its_phase(self):
        self.complete(2)
        self.c.advance("pilot", [["200_3", "RUNNING", "shared", "48:00:00"]],
                       {"100_1": "OUT_OF_MEMORY"})
        self.c.retry.assert_called_once()
        phase, failures, attempt = self.c.retry.call_args[0]
        self.assertEqual((phase, failures, attempt["job"]), ("pilot", {1: "OUT_OF_MEMORY"}, "100"))

    def test_retry_records_reservation_resources_and_phase(self):
        self.c.sbatch = Mock(return_value="400")
        attempt = self.c.state["pilot"]["arrays"][0]
        module.Controller.retry(self.c, "pilot", {1: "OUT_OF_MEMORY"}, attempt)
        args = self.c.sbatch.call_args[0]
        self.assertIn("--array=1%2", args[0])
        self.assertEqual(args[-2:], (48, 32768))
        self.assertEqual((self.study / "submissions.tsv").read_text(), "400\tpilot\t32768\t1\t48\n")
        self.assertEqual(self.c.state["pilot"]["arrays"][-1]["ids"], [1])
        self.assertIsNone(self.c.state["submission_in_progress"])
        with self.assertRaises(ValueError):
            module.Controller.retry(self.c, "pilot", {1: "OUT_OF_MEMORY"}, attempt)

    def test_failed_submission_keeps_reservation(self):
        self.c.sbatch = Mock(side_effect=module.subprocess.CalledProcessError(1, ["sbatch"]))
        with self.assertRaises(module.subprocess.CalledProcessError):
            module.Controller.retry(self.c, "full", {3: "NODE_FAIL"}, self.c.state["full"]["arrays"][0])
        self.assertIn("3", self.c.state["retried"])
        self.assertIsNotNone(self.c.state["submission_in_progress"])
        self.assertEqual(len(self.c.state["full"]["arrays"]), 1)

    def test_validation_submission_is_only_a_validator(self):
        self.c.sbatch = Mock(return_value="500")
        module.Controller.validate(self.c, "pilot")
        args = self.c.sbatch.call_args[0]
        self.assertEqual(args[1].name, "validate-study-stage.sh")
        self.assertEqual(args[2], [str(self.study), "pilot"])
        self.assertEqual((self.study / "controllers.tsv").read_text(), "500\tpilot\t100\n")
        self.assertEqual(self.c.state["pilot"]["validator"], "500")
        self.assertEqual(self.c.sbatch.call_args[1]["hours"], 12)

    def test_full_collection_has_time_for_all_replicates(self):
        self.c.sbatch = Mock(return_value="501")
        module.Controller.validate(self.c, "full")
        self.assertEqual(self.c.sbatch.call_args[1]["hours"], 48)

    def test_routing_prioritizes_validator_within_available_capacity(self):
        self.c.state["pilot"]["validator"] = "300"
        self.c.wall_hours = lambda value: 48
        rows = [["100_1", "RUNNING", "highmem", "48:00:00"],
                ["100_2", "RUNNING", "highmem", "48:00:00"],
                ["200_3", "PENDING", "shared", "48:00:00"],
                ["300", "PENDING", "shared", "12:00:00"],
                ["999_1", "PENDING", "shared", "48:00:00"]]
        module.Controller.route(self.c, rows)
        self.assertEqual([call[0][0][2] for call in self.c.command.call_args_list],
                         ["JobId=300", "JobId=200_3"])

    def test_fatal_failure_is_detected_before_other_tasks_finish(self):
        with self.assertRaises(ValueError):
            self.c.advance("full", [["200_4", "RUNNING", "shared", "48:00:00"]], {"200_3": "FAILED"})
        self.c.retry.assert_not_called()

    def test_validator_exit_alone_cannot_create_acceptance(self):
        self.c.state["pilot"]["validator"] = "300"
        with self.assertRaises(ValueError):
            self.c.advance("pilot", [], {"300": "COMPLETED"})
        self.assertFalse(self.c.state["pilot"]["accepted"])
        (self.study / "PILOT-ACCEPTED").write_text("fixture\n")
        self.c.advance("pilot", [], {"300": "COMPLETED"})
        self.assertTrue(self.c.state["pilot"]["accepted"])

    def test_collector_requires_completed_results_and_pilot_acceptance(self):
        self.complete(3, 4)
        self.c.advance("full", [], {})
        self.c.validate.assert_not_called()
        self.c.state["pilot"]["accepted"] = True
        self.c.advance("full", [], {})
        self.c.validate.assert_called_once_with("full")

    def test_lost_submission_response_blocks_duplicate_retry(self):
        self.c.state["submission_in_progress"] = {"phase": "full", "tasks": [3]}
        self.c.step()
        self.c.block.assert_called_once()
        self.c.retry.assert_not_called()

    def test_throttle_increases_only_after_validated_pilot(self):
        self.complete(1, 2, 3, 4)
        self.c.queue = Mock(return_value=[["200_4", "RUNNING", "shared", "48:00:00"]])
        self.c.accounting = Mock(return_value={})
        self.c.step()
        self.c.command.assert_not_called()
        self.c.state["pilot"]["accepted"] = True
        self.c.step()
        self.c.command.assert_called_once_with(["scontrol", "update", "JobId=200", "ArrayTaskThrottle=32"])


if __name__ == "__main__":
    unittest.main()
