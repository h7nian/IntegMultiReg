"""Verify the archived release's specific leak, regardless of leak classification."""
import pathlib
import re
import sys

text = pathlib.Path(sys.argv[1]).read_text()
assert '[ FAIL 0 | WARN 0 | SKIP 0 | PASS 217 ]' in text, 'Old-release tests did not pass'
pattern = r'(?m)^==\d+== ([\d,]+) bytes in ([\d,]+) blocks are (definitely|possibly|indirectly) lost in loss record[^\n]*\n(.*?)(?=^==\d+== \s*$)'
records = list(re.finditer(pattern, text, re.S | re.M))
assert records, 'Missing detailed leak records'
total_bytes = total_blocks = 0
for record in records:
    size, blocks, kind, stack = record.groups()
    assert kind in ('definitely', 'possibly'), 'Unexpected indirect leak'
    for frame in ('dvector (utils.c:756)', 'predict_cv_fold (prediction_cv.c:28)',
                  'main_function_prediction (main_prediction.c:246)'):
        assert frame in stack, 'Unexpected leak allocation stack: ' + frame
    total_bytes += int(size.replace(',', ''))
    total_blocks += int(blocks.replace(',', ''))
assert (total_bytes, total_blocks) == (7200, 60), (total_bytes, total_blocks)
summary_bytes = summary_blocks = 0
for kind in ('definitely', 'indirectly', 'possibly'):
    matches = re.findall(rf'{kind} lost:\s*([\d,]+) bytes in ([\d,]+) blocks', text)
    assert len(matches) == 1, 'Missing/duplicate leak summary: ' + kind
    size, blocks = (int(x.replace(',', '')) for x in matches[0])
    if kind == 'indirectly':
        assert (size, blocks) == (0, 0), 'Unexpected indirect leak summary'
    summary_bytes += size
    summary_blocks += blocks
assert (summary_bytes, summary_blocks) == (7200, 60), 'Leak records/summary mismatch'
print('Old-release tests passed; verified 7,200 bytes / 60 blocks from the historical prediction allocation stack (definite + possible).')
