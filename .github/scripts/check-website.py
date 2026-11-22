"""Check generated local navigation, reference coverage and required guides."""
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit
import json
import re
import sys


class Page(HTMLParser):
    def __init__(self, text):
        super().__init__(convert_charrefs=True)
        self.ids = set()
        self.links = []
        self.viewport = False
        self.redirect = False
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if attrs.get("id"):
            self.ids.add(attrs["id"])
        if tag == "a" and attrs.get("name"):
            self.ids.add(attrs["name"])
        for key in ("href", "src"):
            if attrs.get(key):
                self.links.append(attrs[key])
        if tag == "meta" and attrs.get("name") == "viewport":
            self.viewport = True
        if tag == "meta" and attrs.get("http-equiv", "").lower() == "refresh":
            match = re.search(r"url\s*=\s*(.+)$", attrs.get("content", ""), re.I)
            if match:
                self.redirect = True
                self.links.append(match.group(1).strip(" \"'"))


root = Path(sys.argv[1]).resolve()
pages = {p: Page(p.read_text(encoding="utf-8")) for p in root.rglob("*.html")}
if not pages:
    raise SystemExit("No website pages were built.")
errors = []
checked = 0
for path, page in pages.items():
    if not page.viewport and not page.redirect:
        errors.append("Missing viewport: " + str(path.relative_to(root)))
    for link in page.links:
        url = urlsplit(link)
        if url.scheme not in ("", "http", "https"):
            continue
        if url.netloc and url.netloc != "h7nian.github.io":
            continue
        target = unquote(url.path)
        if url.netloc or target.startswith("/"):
            if not target.startswith("/IntegMultiReg/"):
                continue
            dest = root / target[len("/IntegMultiReg/"):]
        else:
            dest = path.parent / target if target else path
        dest = dest.resolve()
        if dest.is_dir():
            dest = dest / "index.html"
        checked += 1
        if not dest.exists():
            errors.append(str(path.relative_to(root)) + " -> missing " + link)
        elif url.fragment and dest in pages:
            anchor = unquote(url.fragment)
            if anchor not in pages[dest].ids:
                errors.append(str(path.relative_to(root)) + " -> missing anchor " + link)

required = [
    "index.html", "reference/index.html", "articles/IntegMultiReg.html",
    "articles/uncertainty.html", "articles/interpretation.html",
    "news/index.html", "authors.html",
]
for page in required:
    if not (root / page).is_file():
        errors.append("Missing page: " + page)
package = Path("man")
if package.is_dir():
    for rd in package.glob("*.Rd"):
        if not (root / "reference" / (rd.stem + ".html")).is_file():
            errors.append("Unpublished help topic: " + rd.stem)
home = (root / "index.html").read_text(encoding="utf-8")
if not re.search(r"\b0\.1\.3\b", home):
    errors.append("The homepage does not identify version 0.1.3.")
result = {
    "pages": len(pages),
    "content_pages": sum(not page.redirect for page in pages.values()),
    "redirect_pages": sum(page.redirect for page in pages.values()),
    "local_links_checked": checked,
    "errors": errors,
}
print(json.dumps(result, indent=2))
if errors:
    raise SystemExit(1)
