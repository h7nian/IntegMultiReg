"""Write a read-only progress report and score previews from completed tasks.

Reports never create acceptance markers or combine audit scores with scientific
scores. Large fit/checkpoint RDS files are not loaded by this metadata monitor.
"""
import argparse
from collections import Counter
import csv
import datetime
import fcntl
import io
import json
import os
from pathlib import Path
import re
import subprocess


def atomic_text(path, text):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text)
    temporary.replace(path)


def csv_text(rows, fields):
    out = io.StringIO()
    writer = csv.DictWriter(out, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
    return out.getvalue()


def ledger(path, columns):
    if not path.exists():
        return []
    return [row for row in csv.reader(path.read_text().splitlines(), delimiter="\t")
            if len(row) == columns and re.fullmatch(r"[1-9][0-9]*", row[0])]


def report(study):
    tasks = list(csv.DictReader((study / "tasks.csv").open()))
    submissions = ledger(study / "submissions.tsv", 5)
    controllers = ledger(study / "controllers.tsv", 3)
    latest_jobs, latest_phases = {}, {}
    for job, phase, memory, ids, hours in submissions:
        for task_id in ids.split(","):
            latest_jobs[int(task_id)] = job + "_" + str(int(task_id))
            latest_phases[int(task_id)] = phase
    known = {row[0] for row in submissions + controllers}
    queue_output = subprocess.check_output([
        "squeue", "-h", "-r", "-u", os.environ["USER"],
        "-o", "%i|%T|%P|%M|%R"], universal_newlines=True, timeout=30)
    queue = {}
    for line in queue_output.splitlines():
        fields = [field.strip() for field in line.split("|", 4)]
        if len(fields) == 5 and fields[0].split("_")[0] in known:
            queue[fields[0]] = dict(zip(("job", "state", "partition", "elapsed", "reason"), fields))
    accounting = {}
    if known:
        start = datetime.datetime.fromtimestamp((study / "study.rds").stat().st_mtime - 86400)
        output = subprocess.check_output([
            "sacct", "-nP", "-j", ",".join(sorted(known)),
            "--starttime=" + start.strftime("%Y-%m-%d"),
            "--format=JobID%40,State%40,ExitCode,Elapsed,MaxRSS"],
            universal_newlines=True, timeout=30)
        for line in output.splitlines():
            fields = [field.strip() for field in line.split("|")]
            if len(fields) >= 5 and "." not in fields[0]:
                accounting[fields[0]] = dict(zip(("job", "state", "exit_code", "elapsed", "max_rss"), fields))
    running = {"RUNNING", "COMPLETING", "CONFIGURING", "SUSPENDED", "RESIZING", "STAGE_OUT"}
    failures = {"FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED", "BOOT_FAIL", "DEADLINE"}
    task_rows, rank_audits = [], []
    score_rows = {"audit": [], "scientific": []}
    score_fields = {"audit": set(), "scientific": set()}
    groups = {"audit": Counter(), "pilot": Counter(), "full": Counter()}
    for task in tasks:
        task_id = int(task["task_id"])
        directory = study / "tasks" / ("%04d" % task_id)
        group = "audit" if task["kind"] == "audit" else "pilot" if task["pilot"] == "TRUE" else "full"
        job = latest_jobs.get(task_id, "")
        observed = queue.get(job, accounting.get(job, {}))
        state = observed.get("state", "SUBMITTED" if job else "UNSUBMITTED").split()[0].rstrip("+")
        completed = (directory / "TASK-COMPLETED").exists()
        if completed:
            state, category = "COMPLETED", "completed"
        elif state in running:
            category = "running"
        elif state in failures or state == "COMPLETED":
            category = "attention"
        elif state == "UNSUBMITTED":
            category = "unsubmitted"
        else:
            category = "pending"
        groups[group][category] += 1
        task_rows.append(dict(task, stage=group, job=job, state=state, category=category,
                              submission_phase=latest_phases.get(task_id, ""),
                              partition=observed.get("partition", ""), reason=observed.get("reason", ""),
                              exit_code=observed.get("exit_code", "")))
        if not completed:
            continue
        for path in sorted(directory.glob("configuration-*/*.rank-audit.csv")):
            rank_audits.append({"task_id": task_id, "path": str(path.relative_to(study))})
        bucket = "audit" if group == "audit" else "scientific"
        for path in sorted(directory.glob("configuration-*/*-scores.csv")):
            for score in csv.DictReader(path.open()):
                fields = {"score_" + key: value for key, value in score.items()}
                score_fields[bucket].update(fields)
                score_rows[bucket].append(dict(task, stage=group,
                    method=path.name[:-len("-scores.csv")], **fields))
    current_array = submissions[-1][0] if submissions else ""
    controller_problems = []
    overlap_policy = study / "execution/overlap-policy.json"
    policy = json.loads(overlap_policy.read_text()) if overlap_policy.exists() else {}
    superseded = set(policy.get("superseded_controllers", []))
    last_controller = {phase: (job, array) for job, phase, array in controllers if job not in superseded}
    for phase, (job, array) in last_controller.items():
        state = queue.get(job, accounting.get(job, {})).get("state", "UNKNOWN")
        if state.split()[0].rstrip("+") in failures:
            controller_problems.append({"job": job, "phase": phase, "state": state})
    execution_path = study / "execution/overlap-state.json"
    execution = json.loads(execution_path.read_text()) if execution_path.exists() else {}
    now = datetime.datetime.now(datetime.timezone.utc).astimezone().isoformat(timespec="seconds")
    snapshot = {"updated": now, "current_array": current_array,
        "current_phase": submissions[-1][1] if submissions else "not submitted",
        "groups": {group: {key: counts[key] for key in
                   ("completed", "running", "pending", "unsubmitted", "attention")}
                   for group, counts in groups.items()},
        "queue": list(queue.values()), "rank_audits": rank_audits,
        "controller_problems": controller_problems,
        "execution_blocker": execution.get("blocked"),
        "concurrency_limit": policy.get("concurrency", 16),
        "gates": {name: (study / name).exists() for name in
                  ("AUDITS-ACCEPTED", "PILOT-ACCEPTED", "VALIDATION-COMPLETED")}}
    progress = study / "progress"
    progress.mkdir(exist_ok=True)
    atomic_text(progress / "snapshot.json", json.dumps(snapshot, indent=2) + "\n")
    atomic_text(progress / "tasks.csv", csv_text(task_rows, list(task_rows[0])))
    base_fields = list(tasks[0]) + ["stage", "method"]
    for bucket in score_rows:
        fields = base_fields + sorted(score_fields[bucket])
        atomic_text(progress / (bucket + "-scores.csv"), csv_text(score_rows[bucket], fields))
    names = {"audit": "原法审计", "pilot": "完整参数试跑", "full": "其余全量任务"}
    lines = ["# IntegMultiReg 运行进度", "", "更新时间：" + now, "",
        "最近提交阶段：`" + snapshot["current_phase"] + "`；最近数组：`" + current_array + "`。", "",
        "研究总并发上限：%d；试跑与全量可同时运行。" % snapshot["concurrency_limit"] if policy else
        "研究总并发上限：16。", "",
        "| 任务 | 总数 | 完成 | 运行 | 排队/等待记账 | 未提交 | 待检查 |",
        "|---|---:|---:|---:|---:|---:|---:|"]
    for group, label in names.items():
        counts = groups[group]
        lines.append("| " + label + " | " + " | ".join(str(n) for n in
            [sum(counts.values())] + [counts[key] for key in
             ("completed", "running", "pending", "unsubmitted", "attention")]) + " |")
    lines += ["", "已保存 %d 项零 ridge 秩亏审计；这些记录与完整参数实验分开。" % len(rank_audits), "",
        "- [已完成科学任务的分数](progress/scientific-scores.csv)",
        "- [审计任务分数](progress/audit-scores.csv)",
        "- [逐任务状态](progress/tasks.csv)",
        "- [当前队列和验收标记](progress/snapshot.json)", "",
        "- [软件验证与研究计划](../../IMPLEMENTATION.md)", "",
        "分数只纳入已有 TASK-COMPLETED 的任务，各行保留 reference、ridge 和方法标签。",
        "这些是阶段结果；全研究比较、收敛和稿件结论仍需后续验收。", ""]
    if controller_problems:
        lines += ["当前阶段验收任务需要检查：", ""]
        lines += ["- `{job}`：{phase}，{state}；[日志](logs/accept-{phase}-{job}.log)。".format(**row)
                  for row in controller_problems]
        lines.append("")
    if snapshot["execution_blocker"]:
        lines += ["执行需要检查：" + snapshot["execution_blocker"], "",
                  "详见 [执行记录](execution/ATTENTION.md)。", ""]
    atomic_text(study / "PROGRESS.md", "\n".join(lines))
    print("Progress snapshot refreshed:", now, snapshot["current_phase"], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path)
    args = parser.parse_args()
    with (args.study / ".progress.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        report(args.study)


if __name__ == "__main__":
    main()
