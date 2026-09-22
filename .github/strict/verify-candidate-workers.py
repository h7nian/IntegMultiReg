import pathlib, re, sys
text = pathlib.Path(sys.argv[1]).read_text()
counts = re.findall(r'ERROR SUMMARY:\s*([0-9,]+) errors', text)
assert counts and all(int(x.replace(',', '')) == 0 for x in counts), 'Missing or nonzero parent Memcheck summary'
assert re.search(r'Verified [1-9][0-9]* instrumented worker processes', text), 'Missing worker verification'
print('Instrumented parent and worker suite passed.')
