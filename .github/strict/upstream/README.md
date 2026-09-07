# Upstream graphics-library Valgrind rules

These files are copied without changes from GNOME upstream, not generated
from IntegMultiReg failures. No IntegMultiReg function is suppressed.
Original unfiltered graphics controls are retained in the audit evidence.

Source URLs:
- https://github.com/GNOME/pango/blob/main/pango.supp
- https://github.com/GNOME/glib/blob/main/tools/glib.supp

Explanation of Fontconfig encoded-pointer false positives:
https://gnome.pages.gitlab.gnome.org/librsvg/devel-docs/memory_leaks.html

Content hashes (snapshot retrieved 2026-09-07 UTC):
- pango.supp: 2f93fad0bd8ae72fdcfd3c73cc960c1163ada6475340deadb5c009beb8d3e934
- glib.supp: 26a27e2737a96289356cd4d04291d4644ebf868206afcf4671cd80abfbcd837c

Use Valgrind's default error leak kinds (definite, possible); retain indirect
loss reports as well. Indirect losses are consequences of direct losses,
including the documented encoded-pointer false positives. The exact old
package leak and the independent R-heap control must still be detected.
A clean result means zero *unsuppressed* errors, not that no system-library
allocations were excluded.
