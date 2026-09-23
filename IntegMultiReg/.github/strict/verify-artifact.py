import hashlib, pathlib, subprocess, tarfile
root = pathlib.Path('.')
fixtures = root / '.github/strict/fixtures'
expected = {
 'IntegMultiReg_0.1.1.tar.gz': '6e5c82096a9d0cc62161a810683ccbb6f3c5a748bf89b13432b1ec8976d07745',
 'IntegMultiReg_0.1.3.tar.gz': '46212e3a011986bf9b0e4cc81e75af137f95db3e731218558d1b9a5e6bca3c3b'
}
for name, digest in expected.items():
    assert hashlib.sha256((fixtures/name).read_bytes()).hexdigest() == digest, name
    print(digest, name)
# Refuse a green result for a stale tarball after package source changes.
tracked = subprocess.check_output(['git','ls-files'],text=True).splitlines()
with tarfile.open(fixtures/'IntegMultiReg_0.1.3.tar.gz') as archive:
    members = set(archive.getnames())
    for name in tracked:
        member = 'IntegMultiReg/'+name
        if name == 'DESCRIPTION':
            continue  # R build adds Packaged/Author/Maintainer metadata.
        if member in members:
            packaged = archive.extractfile(member).read()
            checkout = (root/name).read_bytes()
            if name.endswith('.win'):
                packaged = packaged.replace(b'\r\n', b'\n').rstrip(b'\n')
                checkout = checkout.replace(b'\r\n', b'\n').rstrip(b'\n')
            assert packaged == checkout, name
        elif name.startswith(('R/','man/','tests/','src/','data/')):
            raise AssertionError('Missing source file: '+name)
    assert b'Version: 0.1.3' in archive.extractfile('IntegMultiReg/DESCRIPTION').read()
assert '\nVersion: 0.1.3\n' in (root/'DESCRIPTION').read_text()
print('Pinned release source matches checkout.')
