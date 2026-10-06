"""Check generated local links and assets without making network requests."""
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit
import sys


class Page(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self.ids = set()

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if "id" in attrs:
            self.ids.add(attrs["id"])
        if tag == "a" and "name" in attrs:
            self.ids.add(attrs["name"])
        for name in ("href", "src"):
            if name in attrs:
                self.links.append(attrs[name])


root = Path(sys.argv[1] if len(sys.argv) > 1 else "docs").resolve()
pages = {}
for path in root.rglob("*.html"):
    page = Page()
    page.feed(path.read_text(encoding="utf-8"))
    pages[path] = page
failures = []
for source, page in pages.items():
    for link in page.links:
        url = urlsplit(link)
        if url.scheme == "doi":
            failures.append((source.relative_to(root), link + " (use https://doi.org/)"))
        if url.scheme or url.netloc or not url.path:
            continue
        path = unquote(url.path)
        if path.startswith("/IntegMultiReg/"):
            target = root / path[len("/IntegMultiReg/"):]
        elif path.startswith("/"):
            continue
        else:
            target = source.parent / path
        target = target.resolve()
        if target.is_dir():
            target = target / "index.html"
        if not target.exists():
            failures.append((source.relative_to(root), link))
        elif url.fragment and target in pages:
            if unquote(url.fragment) not in pages[target].ids:
                failures.append((source.relative_to(root), link))
if not pages:
    failures.append((root, "No generated HTML pages"))
for page_name, source_path in {"method-coverage.html": "inst/METHOD-COVERAGE.md",
                               "migration.html": "inst/MIGRATION.md"}.items():
    page = pages.get(root / page_name)
    if page is not None and not any("/blob/" in link and link.endswith("/" + source_path)
                                    for link in page.links):
        failures.append((page_name, "Source link must point to " + source_path))
for source, link in sorted(set(failures), key=str):
    print("Broken local link:", source, "->", link)
if failures:
    sys.exit(1)
print("Checked local links, assets and DOI URL schemes in", len(pages), "HTML pages.")
