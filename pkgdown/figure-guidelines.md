# Figure style

The manuscript, response and introductory guide use the same palette and
typographic scale. These settings affect presentation only. Package plotting
defaults and statistical routines remain those of release 0.1.3.

## Sources

- [JSS style guide](https://www.jstatsoft.org/style): annotations should be
  about the size of the caption, or slightly smaller; use sentence case.
- [Nature research figure guide](https://research-figure-guide.nature.com/figures/building-and-exporting-figure-panels/):
  use standard fonts, embed editable text, align panels, avoid excess white
  space and identify series with labels as well as colour. Its 5–7 pt sizing
  is specific to Nature; the JSS article uses the larger caption-based scale.
- [UMN colours](https://brand.umn.edu/colors): current maroon `#6D001F`,
  gold `#FFCC33`, neutral grey `#5A5A5A`. Accessed November 22, 2026.

## Settings used here

- White background and a standard sans serif font. Paper PDFs embed Nimbus
  Sans, a Helvetica-compatible font. Website figures use the device's sans font.
- At the JSS text width of 6.125 inches, most axis labels are about 8–9 pt,
  axis titles 9–10 pt and panel titles 10 pt. Device width, font size and
  stroke width scale together for the three-panel heatmap.
- Trace strokes are about 0.6 pt at final size; axis strokes about 0.6 pt.
  Reference lines use dark grey dashes, which remain visible on white and in
  greyscale. Gold is used for platform bars with dark annotations.
- Retain complete trajectories and the full 0–1 mPIP scale. Heatmaps use a
  light-to-maroon sequential palette. Platform bars use maroon, gold and grey
  in input order; feature labels also name the platform.
- Use shared scale limits across comparable panels. Keep labels out of the
  data area. Omit a redundant platform legend where labels identify platforms.
  Subgroup annotations on gold use dark text.
- Export paper figures directly to vector PDF with embedded fonts. Inspect
  the compiled article and response at their final sizes, including greyscale.

The export helper is distributed with the manuscript replication materials
as `figure-style.R`; guide examples pass styling arguments to public functions.
