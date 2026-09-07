import pathlib, re, sys
base = pathlib.Path(sys.argv[1])
check = base/'candidate/IntegMultiReg.Rcheck'
log = (check/'00check.log').read_text()
assert not re.search(r'\b(?:ERROR|WARNING)\b',log), 'Package check has ERROR/WARNING'
assert (base/'candidate-exit-status.txt').read_text().strip() == '0'
for label in ('checking examples','checking tests','checking re-building of vignette outputs'):
    blocks = re.split(r'(?m)^\* ', log)
    matches = [block for block in blocks if block.startswith(label)]
    assert len(matches) == 1, 'Missing or duplicate stage: '+label
    header = matches[0].splitlines()[0]
    timing_note = (label == 'checking examples' and header.endswith(' NOTE')
                   and 'Examples with CPU (user + system) or elapsed time > 5s' in matches[0])
    assert re.search(r'\bOK\s*$', header) or timing_note, 'Stage did not finish successfully: '+label
summaries=[]
for path in check.rglob('*'):
    if path.is_file() and path.suffix in ('.Rout','.fail','.log'):
        text=path.read_text(errors='replace')
        for count in re.findall(r'ERROR SUMMARY:\s*([0-9,]+) errors',text):
            assert int(count.replace(',','')) == 0, str(path)
            summaries.append(str(path.relative_to(check)))
assert any(path.endswith('-Ex.Rout') for path in summaries), 'Missing example Valgrind summary'
assert any(path.startswith('tests/') and path.endswith('.Rout') for path in summaries), 'Missing test Valgrind summary'
assert 'build_vignettes.log' in summaries, 'Missing vignette rebuild Valgrind summary'
(base/'candidate-valgrind-files.txt').write_text('\n'.join(summaries)+'\n')
print('Candidate package check and Valgrind summaries passed; inspect any NOTEs separately.')
