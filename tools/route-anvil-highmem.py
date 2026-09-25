"""Add highmem eligibility to one study array within Anvil's four-job limit.

This changes only pending task partition eligibility. Slurm still enforces
the two-running-job highmem limit and the array's existing concurrency cap.
"""
import argparse
import datetime
import fcntl
import json
import os
from pathlib import Path
import re
import subprocess
import time


def wall_hours(value):
    if not re.fullmatch(r"(?:[0-9]+-)?[0-9:]+", value):
        return float("inf")
    days, clock = value.split("-", 1) if "-" in value else ("0", value)
    fields = [int(part) for part in clock.split(":")]
    if len(fields) > 3:
        return float("inf")
    hours, minutes, seconds = [0] * (3 - len(fields)) + fields
    return int(days) * 24 + hours + minutes / 60 + seconds / 3600


def routing_plan(rows, array):
    occupied = sum("highmem" in row[2].split(",") for row in rows)
    waiting = [row for row in rows
               if re.fullmatch(re.escape(array) + r"_[0-9]+", row[0])
               and row[1] == "PENDING" and row[2] == "shared"
               and wall_hours(row[3]) <= 48]
    waiting.sort(key=lambda row: int(row[0].split("_")[1]))
    return waiting, waiting[:max(0, 4 - occupied)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path)
    parser.add_argument("array")
    parser.add_argument("--once", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if not re.fullmatch(r"[1-9][0-9]*", args.array):
        parser.error("array must be a Slurm job ID")
    submissions = (args.study / "submissions.tsv").read_text().splitlines()
    if args.array not in {line.split("\t")[0] for line in submissions}:
        parser.error("array is not recorded in this study")
    with (args.study / (".highmem-routing-" + args.array + ".lock")).open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        print("Routing array", args.array, "within four highmem submissions; wall limit 48 hours.", flush=True)
        while True:
            output = subprocess.check_output([
                "squeue", "-h", "-r", "-u", os.environ["USER"],
                "-o", "%i|%T|%P|%l"], universal_newlines=True)
            rows = [[field.strip() for field in line.split("|")]
                    for line in output.splitlines() if line.strip()]
            assert all(len(row) == 4 for row in rows)
            waiting, requests = routing_plan(rows, args.array)
            for row in requests:
                record = {"time_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                          "job": row[0], "partitions": ["shared", "highmem"]}
                if args.dry_run:
                    record["dry_run"] = True
                else:
                    subprocess.check_call(["scontrol", "update", "JobId=" + row[0],
                                           "Partition=shared,highmem"])
                    with (args.study / "partition-routing-events.jsonl").open("a") as log:
                        log.write(json.dumps(record) + "\n")
                print(json.dumps(record), flush=True)
            if args.once or args.dry_run or not waiting:
                break
            time.sleep(30)


if __name__ == "__main__":
    main()
