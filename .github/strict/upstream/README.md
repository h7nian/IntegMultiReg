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

ubuntu-vignette-font-baseline.supp contains 162 bounded possible-loss patterns
matched to package-free rmarkdown graphics controls from runs 35813024084,
35818036159, 35818897822, 35827733091, 35835326668 and the standalone PNG control in 35837248352, plus scoped controls in 35837968029. The adjacent JSON pins each source-log hash and the
rule hash. Allocation frames and within-object relative call-site addresses
match a control through every retained frame. Most patterns have 24 frames;
ten use 20–23 frames before differing R evaluation context, and three retain
complete 9/10-frame background-thread stacks. No arbitrary-depth or object-name
wildcard is used. Named functions and Ubuntu shared-object paths are exact.

The rules statically cover 772 observed possible-loss reports across controls
and four candidate runs. None matches the archived package leak records.
These are graphics-environment exclusions, not proof of leak-free graphics
libraries or general package correctness. No definite-loss, conditional-jump,
invalid-read/write or general R allocation rule is added here.

Raw controls remain recorded before activation. Standalone PNG, default and width-first
graphics must pass repeated scoped runs, followed by the R-heap read probe and
archived 7,200-byte package-leak negative controls under the same rules. Complete
candidate and worker-suite checks remain mandatory. Raw-control classification
varies between executions, so a single clean replay is insufficient evidence.

Raw standalone PNG and vignette logs are checked before scoped replay: every
error must be a possible-loss context matching this bundle; missing summaries,
other error classes, nonmatching allocations and incorrect exits fail the gate.
Each of standalone PNG, default vignette and width-first vignette then requires
three zero-exit, zero-error scoped replays. Recording an expected raw nonzero
exit does not bypass this gate or the negative controls.
