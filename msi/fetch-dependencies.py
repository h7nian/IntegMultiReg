"""Download exact CRAN source versions from the lock; never substitute versions."""
import csv, hashlib, json, pathlib, sys, urllib.request, time
from concurrent.futures import ThreadPoolExecutor
folder = pathlib.Path(sys.argv[1])
rows = list(csv.DictReader((folder/'dependencies.tsv').open(), delimiter='\t'))
def fetch(row):
    name, version = row['Package'], row['Version']
    target = folder/f'{name}_{version}.tar.gz'
    urls = [f'https://cran.r-project.org/src/contrib/{target.name}',
            f'https://cran.r-project.org/src/contrib/Archive/{name}/{target.name}']
    errors = []
    for url in urls:
        for attempt in range(2):
            try:
                with urllib.request.urlopen(url, timeout=90) as response:
                    data = response.read()
                if data[:2] != b'\x1f\x8b': raise ValueError('Not gzip')
                target.write_bytes(data)
                result = dict(package=name, version=version, file=target.name, url=url,
                              bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
                print(f'{name} {version}: {len(data)} bytes', flush=True)
                return result
            except Exception as exc:
                errors.append(f'{url}: {exc}')
                if attempt == 0: time.sleep(1)
    raise RuntimeError('\n'.join(errors))
with ThreadPoolExecutor(max_workers=4) as pool:
    results = list(pool.map(fetch, rows))
(folder/'sources.json').write_text(json.dumps(results, indent=2)+'\n')
