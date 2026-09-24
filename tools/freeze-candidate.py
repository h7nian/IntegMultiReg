"""Freeze package sources for an isolated validation build; never overwrite."""
from pathlib import Path
import hashlib, shutil, sys
root = Path(__file__).resolve().parents[1]
out = root/'runtime'/sys.argv[1]
if out.exists():
    raise SystemExit('Choose a new candidate identifier')
out.mkdir(parents=True)
shutil.copytree(root/'IntegMultiReg', out/'source',
    ignore=shutil.ignore_patterns('.git', '*.o', '*.so', '*.dll',
                                 'config.log', 'config.status', 'Makevars', 'Makevars 2'))
for name in ('configure', 'cleanup'):
    (out/'source'/name).chmod(0o755)
lines = []
for path in sorted((out/'source').rglob('*')):
    if path.is_file():
        lines.append(hashlib.sha256(path.read_bytes()).hexdigest()+'  '+str(path.relative_to(out)))
(out/'SOURCE-SHA256SUMS').write_text('\n'.join(lines)+'\n')
print(out)
