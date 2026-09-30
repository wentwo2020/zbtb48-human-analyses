# Figure panel rendering

Two offline scripts that redraw published panels and recompute every statistic
shown on them, from deposited source tables and with no network access. They
exist so that a reader can regenerate a panel and check each plotted value,
without re-running the upstream API queries.

Both read their inputs from the working directory rather than from `../results`,
so copy the tables you need next to the script, or run from a directory holding
both.

## `make_fig5c_panel.py`

Redraws Figure 5C and 5D and recomputes the paired statistics.

Input: `fig5c_donor_by_celltype_matrix.csv`, written by
`../src/04_snrna_cellxgene_spinal.py`.

Output: `fig5c_regenerated.png`, `fig5c_statistics_regenerated.csv`.

The statistics are recomputed from the donor matrix rather than read from a
file, so the numbers on the regenerated panel are derived, not copied.

## `regraph_fig5_panels.py`

Redraws the Figure 5H enrichment panel and the two Wald-ratio forest panels of
Extended Data Figure 7.

Input: `panel_h_source_data.csv`, `panel_i_source_data.csv`,
`panel_j_source_data.csv` and `figure_display_labels.csv`. The first three are
subsets of the table written by `../src/06_finngen_mr_and_enrichment.py`; the
label file maps FinnGen endpoint strings to the short axis labels used in the
figure, which are abbreviations chosen for the figure and are not the source
database's own strings.

Output: `fig5_pain_specificity_REGRAPHED.png`.

## On reproducing the raster exactly

The regenerated images are not pixel-identical to the published ones. Both
scripts save with a tight bounding box, which crops to the rendered extent of
the text, and that extent depends on font metrics resolved at render time. Every
plotted coordinate and every axis string matches; the canvas differs by a few
pixels. If a byte-identical raster is needed, replace the tight bounding box
with an explicit figure size and margins.

One trap worth naming: endpoints for which FinnGen releases no effect estimate
must not be drawn at the null value, because a marker at 1.0 reads as a measured
null rather than as a suppressed estimate. They are omitted from the effect axis
and named in the caption instead.
