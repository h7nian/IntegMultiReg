import pathlib, re, sys
base = pathlib.Path(sys.argv[1])
check = base/'candidate/IntegMultiReg.Rcheck'
log = (check/'00check.log').read_text()
assert not re.search(r'\b(?:ERROR|WARNING)\b',log), 'Package check has ERROR/WARNING'
assert (base/'candidate-exit-status.txt').read_text().strip() == '0'
for label in ('checking examples','checking tests','checking re-building of vignette outputs'):
    assert label in log, 'Missing stage: '+label
summaries=[]
for path in check.rglob('*'):
    if path.is_file() and path.suffix in ('.Rout','.fail','.log'):
        text=path.read_text(errors='replace')
        for count in re.findall(r'ERROR SUMMARY:\s*([0-9,]+) errors',text):
            assert int(count.replace(',','')) == 0, str(path)
            summaries.append(str(path.relative_to(check)))
assert summaries, 'No Valgrind summaries found'
(base/'candidate-valgrind-files.txt').write_text('\n'.join(summaries)+'\n')
print('Candidate package check and Valgrind summaries passed; inspect any NOTEs separately.')
