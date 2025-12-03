# Mantel Correlation Overlay

A small R workflow to compute Pearson correlations and Mantel tests between agronomic variables and visualize them as a heatmap overlaid with Mantel edges.

## Files
- `cor_mantel_overlay.R`: Main script. Builds a Pearson upper-triangular heatmap with significance stars, computes Mantel tests for selected targets, overlays Mantel edges, and saves `cor_upper_heatmap_mantel.pdf` and `.svg`.
- `copy_data_analysis.R`: Additional analysis helper (kept as-is).
- `demo.csv`: Input data expected by the scripts.

## Prerequisites
- R >= 4.0
- Packages: `ggplot2`, `vegan`, `ggnewscale`, `svglite`

The script auto-installs missing packages from `cloud.r-project.org` if they are absent.

## How to run
```bash
Rscript cor_mantel_overlay.R
```

Outputs will appear in the repository root as:
- `cor_upper_heatmap_mantel.pdf`
- `cor_upper_heatmap_mantel.svg`

## What the plot shows
- **Heatmap**: Pearson correlation `r` between variables (upper triangle), colored blue-white-red and annotated with significance stars.
- **Significance legend**: `***` for p < 0.001, `**` for p < 0.01, `*` for p < 0.05.
- **Mantel overlays**: Lines from target nodes (`Yield`, `PUE`, `SB`) to variables.
  - Color encodes Mantel p (orange: 0.001–0.01; blue: 0.01–0.05; grey: > 0.05).
  - Linewidth encodes Mantel r magnitude (thicker for larger r).
  - Linetype encodes sign (solid for positive, dashed for negative).

## Repro tips
- Keep the column names in `demo.csv` unchanged; the script selects specific columns.
- If adding new targets or variables, update the `targets` and `vars` vectors in `cor_mantel_overlay.R`.
- For consistent randomization in Mantel permutations, the seed is set to `123`.
