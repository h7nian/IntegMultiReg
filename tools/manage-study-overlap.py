"""Run production alongside an existing pilot, preserving validation gates.

Scientific task sources and the study manifest remain immutable. A reserved
pilot/retry budget and the production throttle share the configured total cap.
Production uses the full cap after pilot acceptance. Only terminal resource failures
receive one bounded retry. Numerical/cancelled failures stop new admissions.
"""
import argparse
import datetime
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time

RESOURCE_FAILURES = {"OUT_OF_MEMORY", "TIMEOUT", "NODE_FAIL", "PREEMPTED"}
TERMINAL = RESOURCE_FAILURES | {"COMPLETED", "FAILED", "CANCELLED", "BOOT_FAIL", "DEADLINE"}


def atomic_json(path, value):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def retry_request(failures, memory, hours, retried):
    if not failures or any(state not in RESOURCE_FAILURES or str(task) in retried
                           for task, state in failures.items()):
        raise ValueError("Unexplained/cancelled failure or exhausted resource retry")
    memory *= 2 if "OUT_OF_MEMORY" in failures.values() else 1
    hours *= 2 if "TIMEOUT" in failures.values() else 1
    if hours > 96 or math.ceil(memory / 1896) > 128:
        raise ValueError("Resource retry exceeds the recorded Anvil limits")
    return memory, hours


class Controller:
    def __init__(self, study):
        self.study = study.resolve()
        self.directory = self.study / "execution"
        self.policy = json.loads((self.directory / "overlap-policy.json").read_text())
        self.state_path = self.directory / "overlap-state.json"
        self.state = json.loads(self.state_path.read_text())
        if not 1 <= self.policy["pilot_reserve"] < self.policy["concurrency"]:
            raise ValueError("Invalid concurrency reservation")
        spec = importlib.util.spec_from_file_location("routing", self.directory / "route-anvil-highmem.py")
        routing = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(routing)
        self.wall_hours = routing.wall_hours
        self.verify()

    def verify(self):
        for name, expected in self.policy["execution_source_hashes"].items():
            if hashlib.sha256((self.directory / name).read_bytes()).hexdigest() != expected:
                raise ValueError("Execution source changed: " + name)
        if hashlib.md5((self.study / "study.rds").read_bytes()).hexdigest() != self.policy["manifest_md5"]:
            raise ValueError("Scientific study manifest changed")
        for name in ("UNIT-VALIDATED", "SANITIZER-VALIDATED", "NATIVE-VALIDATED"):
            if not (Path(self.policy["candidate"]) / name).exists():
                raise ValueError("Missing candidate acceptance: " + name)
        if not (self.study / "AUDITS-ACCEPTED").exists():
            raise ValueError("Strict-method audit acceptance is missing")

    def save(self):
        self.state["updated_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        atomic_json(self.state_path, self.state)

    def event(self, message):
        print(datetime.datetime.now(datetime.timezone.utc).isoformat(), message, flush=True)

    def command(self, arguments):
        return subprocess.check_output(arguments, universal_newlines=True, timeout=45).strip()

    def queue(self):
        output = self.command(["squeue", "-h", "-r", "-u", os.environ["USER"],
                               "-o", "%i|%T|%P|%l"])
        rows = [[field.strip() for field in line.split("|")] for line in output.splitlines() if line]
        if any(len(row) != 4 for row in rows):
            raise ValueError("Unexpected Slurm queue record")
        return rows

    def accounting(self):
        ids = [item["job"] for phase in ("pilot", "full") for item in self.state[phase]["arrays"]]
        ids += [self.state[p]["validator"] for p in ("pilot", "full") if self.state[p]["validator"]]
        output = self.command(["sacct", "-nP", "-j", ",".join(ids),
                               "--starttime=" + self.policy["accounting_start"],
                               "--format=JobID%40,State%40,ExitCode"])
        states = {}
        for line in output.splitlines():
            fields = line.split("|")
            if len(fields) >= 3 and "." not in fields[0] and fields[1].strip():
                state = fields[1].strip().split()[0].rstrip("+")
                if state == "COMPLETED" and fields[2].strip() != "0:0":
                    state = "FAILED"
                states[fields[0].strip()] = state
        return states

    def block(self, reason):
        self.state["blocked"] = reason
        self.save()
        holds = {}
        for phase in ("pilot", "full"):
            job = self.state[phase]["arrays"][-1]["job"]
            try:
                result = subprocess.run(["scontrol", "hold", job], stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE, timeout=30)
                holds[job] = result.returncode
            except subprocess.TimeoutExpired:
                holds[job] = "timeout"
        self.state["hold_requests"] = holds
        self.save()
        (self.directory / "ATTENTION.md").write_text(
            "# Execution needs review\n\n" + reason +
            "\n\nHolds were requested for pending tasks; inspect hold_requests in overlap-state.json."
            " Running tasks and results were retained.\n")
        self.event("BLOCKED: " + reason)

    def sbatch(self, extra, script, arguments, name, hours=12, memory=16384):
        cpus = int(math.ceil(memory / 1896))
        command = ["sbatch", "--parsable", "--account=" + self.policy["account"],
            "--partition=shared", "--nodes=1", "--ntasks=1", "--cpus-per-task=" + str(cpus),
            "--mem=" + str(memory) + "M", "--time=" + str(hours) + ":00:00",
            "--job-name=" + name, "--chdir=" + str(self.study / "source")]
        job = self.command(command + extra + [str(script)] + arguments).split(";")[0]
        if not re.fullmatch(r"[1-9][0-9]*", job):
            raise ValueError("Unexpected sbatch response")
        return job

    def retry(self, phase, failures, attempt):
        memory, hours = retry_request(failures, attempt["memory_mib"], attempt["hours"], self.state["retried"])
        ids = sorted(failures)
        if phase == "pilot" and len(ids) > self.policy["pilot_reserve"]:
            raise ValueError("Pilot retries exceed their concurrency reservation")
        for task in ids:
            self.state["retried"][str(task)] = {"phase": phase, "state": failures[task]}
        # Reserve before submission, so a lost response cannot duplicate a retry.
        self.state["submission_in_progress"] = {"phase": phase, "tasks": ids}
        self.save()
        slots = self.policy["pilot_reserve"] if phase == "pilot" else self.production_slots()
        job = self.sbatch(["--array=" + ",".join(map(str, ids)) + "%" + str(slots),
            "--output=" + str(self.study / "logs/%x-%A_%a.log")],
            self.study / "source/hpc/task.sh", [str(self.study)], "imr-" + phase, hours, memory)
        self.state["submission_in_progress"]["job"] = job
        self.save()
        with (self.study / "submissions.tsv").open("a") as stream:
            stream.write("\t".join([job, phase, str(memory), ",".join(map(str, ids)), str(hours)]) + "\n")
        (self.study / "last-array-id").write_text(job + "\n")
        self.state[phase]["arrays"].append(dict(job=job, ids=ids, memory_mib=memory, hours=hours))
        self.state["submission_in_progress"] = None
        self.save()
        self.event("Submitted bounded " + phase + " retry " + job)

    def validate(self, phase):
        self.state["submission_in_progress"] = {"phase": phase, "validation": True}
        self.save()
        job = self.sbatch(["--output=" + str(self.study / ("logs/accept-" + phase + "-%j.log"))],
            self.directory / "validate-study-stage.sh", [str(self.study), phase], "imr-accept-" + phase)
        self.state["submission_in_progress"]["job"] = job
        self.save()
        array = self.state[phase]["arrays"][-1]["job"]
        with (self.study / "controllers.tsv").open("a") as stream:
            stream.write("\t".join([job, phase, array]) + "\n")
        self.state[phase]["validator"] = job
        self.state["submission_in_progress"] = None
        self.save()
        self.event("Submitted " + phase + " validation " + job)

    def production_slots(self):
        return self.policy["concurrency"] - (0 if self.state["pilot"]["accepted"]
                                             else self.policy["pilot_reserve"])

    def advance(self, phase, rows, states):
        current = self.state[phase]
        if current["accepted"]:
            return
        if current["validator"]:
            state = states.get(current["validator"])
            if state == "COMPLETED":
                marker = "PILOT-ACCEPTED" if phase == "pilot" else "VALIDATION-COMPLETED"
                if not (self.study / marker).exists():
                    raise ValueError(phase + " validator exited without its acceptance evidence")
                current["accepted"] = True
                self.save()
                self.event(phase + " validation passed")
            elif state in TERMINAL:
                raise ValueError(phase + " validation failed: " + str(state))
            return
        latest = {}
        for attempt in current["arrays"]:
            latest.update({task: attempt for task in attempt["ids"]})
        missing = [task for task in self.policy[phase + "_ids"]
                   if not (self.study / "tasks" / ("%04d" % task) / "TASK-COMPLETED").exists()]
        for task in missing:
            state = states.get(latest[task]["job"] + "_" + str(task))
            if state in TERMINAL - RESOURCE_FAILURES - {"COMPLETED"}:
                raise ValueError("Task %d stopped with %s; no automatic numerical/cancel retry" % (task, state))
        arrays = {attempt["job"] for attempt in current["arrays"]}
        if any(row[0].split("_")[0] in arrays for row in rows):
            return
        if missing:
            failures = {task: states.get(latest[task]["job"] + "_" + str(task)) for task in missing}
            if any(state not in TERMINAL for state in failures.values()):
                return  # Accounting may lag behind squeue.
            attempts = {latest[task]["job"] for task in missing}
            if len(attempts) != 1:
                raise ValueError("Missing tasks span multiple exhausted attempts")
            self.retry(phase, failures, latest[missing[0]])
        elif phase == "pilot" or self.state["pilot"]["accepted"]:
            self.validate(phase)

    def route(self, rows):
        occupied = sum("highmem" in row[2].split(",") for row in rows)
        priority = []
        for phase in ("pilot", "full"):
            job = self.state[phase]["validator"]
            priority += [row for row in rows if row[0] == job]
        for phase in ("pilot", "full"):
            arrays = {attempt["job"] for attempt in self.state[phase]["arrays"]}
            priority += sorted([row for row in rows if row[0].split("_")[0] in arrays and "_" in row[0]],
                               key=lambda row: int(row[0].split("_")[1]))
        pending = [row for row in priority if row[1] == "PENDING" and row[2] == "shared"
                   and self.wall_hours(row[3]) <= 48]
        for row in pending[:max(0, 4 - occupied)]:
            self.command(["scontrol", "update", "JobId=" + row[0], "Partition=shared,highmem"])
            self.event("Added highmem eligibility: " + row[0])

    def step(self):
        if self.state.get("blocked"):
            return
        if self.state.get("submission_in_progress"):
            self.block("A prior submission has an unresolved response; inspect before resubmitting")
            return
        try:
            self.verify()
            rows, states = self.queue(), self.accounting()
            for phase in ("pilot", "full"):
                self.advance(phase, rows, states)
            slots = self.production_slots()
            if self.state["full"]["throttle"] != slots:
                job = self.state["full"]["arrays"][-1]["job"]
                if any(row[0].split("_")[0] == job for row in rows):
                    self.command(["scontrol", "update", "JobId=" + job, "ArrayTaskThrottle=" + str(slots)])
                self.state["full"]["throttle"] = slots
                self.save()
            self.route(rows)
        except ValueError as error:
            self.block(str(error))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path)
    parser.add_argument("--once", action="store_true")
    args = parser.parse_args()
    controller = Controller(args.study)
    with (args.study / ".highmem-routing.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        last_report = float("-inf")
        while True:
            try:
                controller.step()
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
                controller.event("Slurm request failed; next poll will refresh state: " + str(error))
            if time.monotonic() - last_report >= 900 or controller.state.get("blocked"):
                try:
                    subprocess.check_call([sys.executable, str(controller.directory / "write-study-progress.py"),
                                           str(controller.study)], timeout=60)
                except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
                    controller.event("Progress refresh failed: " + str(error))
                last_report = time.monotonic()
            if args.once or controller.state.get("blocked") or controller.state["full"]["accepted"]:
                break
            time.sleep(60)


if __name__ == "__main__":
    main()
