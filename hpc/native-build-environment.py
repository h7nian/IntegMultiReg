"""Use the development headers matching Anvil's installed R dependencies."""
from pathlib import Path
import os, re, shlex, subprocess
rhome = Path(subprocess.check_output(['R', 'RHOME'], text=True).strip())
libraries = [rhome/'lib/libR.so', rhome/'library/grDevices/libs/cairo.so']
curl_prefix = Path(subprocess.check_output(['curl-config', '--prefix'], text=True).strip())
prefixes = {curl_prefix}
for library in libraries:
    output = subprocess.check_output(['ldd', str(library)], text=True)
    for path in re.findall(r'(/apps/spack/\S+)', output):
        path = Path(path).resolve()
        if '/apps/r/' not in str(path) and path.parent.name in ('lib', 'lib64'):
            prefixes.add(path.parent.parent)
variables = {'CPATH': [], 'LIBRARY_PATH': [], 'LD_LIBRARY_PATH': [], 'PKG_CONFIG_PATH': []}
for prefix in sorted(prefixes):
    if (prefix/'include').is_dir(): variables['CPATH'].append(str(prefix/'include'))
    for name in ('lib', 'lib64'):
        directory = prefix/name
        if directory.is_dir():
            variables['LIBRARY_PATH'].append(str(directory))
            variables['LD_LIBRARY_PATH'].append(str(directory))
            if (directory/'pkgconfig').is_dir():
                variables['PKG_CONFIG_PATH'].append(str(directory/'pkgconfig'))
    if (prefix/'share/pkgconfig').is_dir():
        variables['PKG_CONFIG_PATH'].append(str(prefix/'share/pkgconfig'))
for key, values in variables.items():
    value = ':'.join(values + ([os.environ[key]] if os.environ.get(key) else []))
    print('export '+key+'='+shlex.quote(value))
