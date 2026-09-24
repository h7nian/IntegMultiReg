"""Derive bounded rules only from the nine completed package-free controls."""
from pathlib import Path
import hashlib, json, re, sys
base = Path(sys.argv[1]).resolve()
assert (base/'COLLECTED').is_file()
logs = sorted(base.glob('*.log'))
assert len(logs) == 9
rules = {}
inputs = {}
for path in logs:
    text = path.read_text()
    assert not re.search(r'(?m)^Error(?:[: ])|Execution halted', text), path
    assert int(Path(str(path)+'.exitcode').read_text()) in (0, 97)
    assert len(re.findall(r'stopifnot\(![\"\']IntegMultiReg[\"\']\s*%in%\s*loadedNamespaces\(\)\)', text)) >= 2
    blocks = re.findall(r'(?m)^\{\n.*?^\}', text, re.S)
    summary = re.findall(r'ERROR SUMMARY:\s*([\d,]+) errors from ([\d,]+) contexts', text)
    assert len(summary) == 1 and int(summary[0][1].replace(',', '')) == len(blocks)
    inputs[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    for block in blocks:
        lines = [line.strip() for line in block.splitlines()]
        assert 'Memcheck:Leak' in lines
        kind = next(line for line in lines if line.startswith('match-leak-kinds:'))
        frames = [line for line in lines if line.startswith(('fun:', 'obj:'))]
        assert len(frames) >= 24 and not any('*' in line or '?' in line for line in frames)
        assert '...' not in lines and not any('imr_' in line.lower() or 'integmultireg' in line.lower() for line in frames)
        language_cache = frames[:3] == ['fun:calloc', 'fun:g_malloc0', 'fun:pango_language_from_string']
        font_pattern = frames[:2] == ['fun:realloc', 'fun:FcPatternObjectInsertElt']
        assert language_cache or font_pattern, frames[:5]
        assert kind == ('match-leak-kinds: possible' if language_cache else 'match-leak-kinds: definite')
        key = (kind, tuple(frames[:24]))
        rules.setdefault(key, []).append(path.name)
output = base/'anvil-font-baseline.supp'
assert not output.exists()
blocks = []
manifest = []
for index, ((kind, frames), sources) in enumerate(sorted(rules.items()), 1):
    name = 'anvil_graphics_baseline_%03d' % index
    blocks.append('{\n   '+name+'\n   Memcheck:Leak\n   '+kind+'\n   '+'\n   '.join(frames)+'\n}\n')
    manifest.append(dict(name=name, leak_kind=kind, frames=list(frames), controls=sources))
output.write_text('\n'.join(blocks))
(base/'rules-audit.json').write_text(json.dumps(dict(
    source_logs=inputs, rules=manifest, rule_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
    scope='Local Fontconfig/Pango baseline only. No wildcard/depth matching, invalid-access or uninitialized-read rules. Requires clean control replays and negative controls before candidate acceptance.'), indent=2)+'\n')
print('Prepared',len(rules),'bounded 24-frame rules from',len(logs),'package-free controls; not yet accepted.')
