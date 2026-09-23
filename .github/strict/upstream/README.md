# Upstream graphics-library Valgrind rules

pango.supp and glib.supp are copied without changes from GNOME upstream, not generated
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

## Audited Ubuntu baseline adaptations

ubuntu-font-baseline.supp is locally authored, NOT claimed to be an upstream
file. It covers the ten remaining contexts from a standalone base-R png/plot
control which explicitly asserts IntegMultiReg is not loaded. The unfiltered
control, upstream-only result and this rule file are retained together.
The one definite-loss pattern is the stripped-symbol counterpart of the
upstream FcPatternObjectInsertElt rule (Fontconfig XML configuration cache).
The other rules are limited to possible losses in Fontconfig/Pango font
caches, GLib reference-count boxes and the Pango background thread TLS.
They do not establish that every system-library allocation is a false
positive; they isolate recorded graphics-environment behavior from package
validation. No invalid-read/write or uninitialized-read error is covered by
these local rules. No IntegMultiReg frame or general R allocation is covered.
Any new unmatched error still fails. The old 7,200-byte package leak and
R-heap uninitialized-read control must continue to fail under these rules.

Baseline evidence run: 34078179119, standalone graphics, ten remaining
contexts after upstream rules, 256 definite and 1,664 possible bytes.

## Audited vignette baseline rules

ubuntu-vignette-font-baseline.supp is a locally generated, control-derived set
of 44 possible-loss patterns. Each contains 24 explicit frames, with exact
function names or exact Ubuntu library paths, and no ellipsis/object wildcard.
Every pattern is observed in the package-free rmarkdown graphics control from
run 35813024084. Its 80 corresponding reports in earlier candidate run
35686577583 match control frames and relative shared-object addresses through
these 24 frames. Two traces differ only much deeper in the R evaluator
(frames 43/44), after the common graphics-library allocation context.

The paired JSON records evidence hashes, rule hash, scope and pending status.
The workflow records the package-free control before activating these additional
rules, then requires its scoped rerun to return zero errors. Before acceptance, rerun the
package-free graphics control, the R-heap uninitialized-read probe and archived
package leak control with the rules active, then the exact candidate and worker
checks. No conditional-jump/invalid-read/write suppression is introduced.
Matched possible-loss graphics allocations are an environment exclusion; this
is not a proof of leak-free graphics libraries or general package correctness.
