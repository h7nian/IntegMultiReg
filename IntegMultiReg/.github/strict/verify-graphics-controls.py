"""Require every raw graphics error to match the audited possible-loss bundle."""
import pathlib
import re
import sys


def patterns(path):
    result = []
    for block in re.findall(r"\{\n(.*?)\n\}", path.read_text(), re.S):
        lines = [line.strip() for line in block.splitlines()]
        assert lines[1:3] == ["Memcheck:Leak", "match-leak-kinds: possible"]
        frames = lines[3:]
        assert frames and all(not any(x in frame for x in ("*", "?", "...")) for frame in frames)
        result.append(frames)
    assert result
    return result


def verify(log, status, rules):
    assert status.read_text().strip() in ("0", "97"), f"Unexpected process exit: {status}"
    text = log.read_text()
    summaries = re.findall(r"ERROR SUMMARY: ([\d,]+) errors from ([\d,]+) contexts", text)
    assert len(summaries) == 1, f"Missing/ambiguous Memcheck summary: {log}"
    records = re.findall(
        r"^==\d+== [\d,]+ bytes in [\d,]+ blocks are (definitely|possibly) lost in loss record[^\n]*\n(.*?)(?=^==\d+== \s*$)",
        text, re.M | re.S,
    )
    errors, contexts = (int(x.replace(",", "")) for x in summaries[0])
    assert status.read_text().strip() == ("97" if errors else "0"), f"Error exit mismatch: {log}"
    assert errors == contexts == len(records), f"Unexpected error class/count: {log}"
    for kind, stack in records:
        assert kind == "possibly", f"Unscoped definite loss: {log}"
        frames = []
        for frame in re.findall(r"^==\d+==\s+(?:at|by) 0x[0-9A-F]+: (.*)$", stack, re.M):
            if frame.startswith("??? (in "):
                frames.append("obj:" + re.fullmatch(r"\?\?\? \(in ([^)]+)\)", frame)[1])
            else:
                frames.append("fun:" + frame.split(" (")[0])
        assert any(frames[:len(rule)] == rule for rule in rules), f"Unmatched raw allocation: {log}"
    print(f"PASS raw graphics: {log.name}, {len(records)} audited possible-loss contexts")


if __name__ == "__main__":
    root, bundle = map(pathlib.Path, sys.argv[1:])
    rules = patterns(bundle)
    logs = [("graphics-after", "graphics-after-exit-status")]
    logs += [("graphics-vignette-control", "graphics-vignette-control-status")]
    logs += [(f"graphics-vignette-control-width-{i}", f"graphics-vignette-control-width-{i}-status") for i in (1, 2, 3)]
    for log, status in logs:
        verify(root / (log + ".log"), root / (status + ".txt"), rules)
