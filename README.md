# mantel-cor_net

English | [中文](README.zh-CN.md)

One R script that draws the "correlation heatmap + Mantel links" figure often seen in soil and crop papers. The upper triangle shows Pearson correlations among explanatory variables. Lines from each response variable (yield, phosphorus use efficiency and so on) to every explanatory variable show a Mantel test, with colour for p, width for |r| and dashes for negative r.

I made it for a phosphorus trial where I wanted soil available P by depth, root and plant P traits, yield and PUE in one figure. The usual tool for this (linkET) was not on CRAN at the time and I needed full control over the layout, so it is plain ggplot2.

## Requirements

R 4.0 or newer with `ggplot2`, `vegan`, `ggnewscale` and `svglite`. Missing packages are installed from CRAN on the first run.

## Usage

```bash
Rscript cor_mantel_overlay.R                    # uses demo.csv
Rscript cor_mantel_overlay.R mydata.csv fig3    # writes fig3.pdf, fig3.svg, fig3_mantel.csv
```

Set the variable lists at the top of the script:

```r
vars    <- c("Year", "cultivar", "Treatment", "0-10cm AP", ...)   # heatmap and link targets
targets <- c("Yield", "PUE", "SB")                                 # response variables
```

`demo.csv` holds 288 plots × 22 columns in the layout the script expects. The values are made up for the demo.

## What is computed

- **Heatmap**: Pearson r and `cor.test` p for each pair of numeric variables, with `*` p < 0.05, `**` p < 0.01, `***` p < 0.001. Pairs with a categorical variable (cultivar, treatment) stay empty, because r on factor codes has no meaning.
- **Mantel links**: for each target and each variable, `vegan::mantel` (Pearson, 999 permutations, seed 123) between two distance matrices. Numeric variables use |a − b|; categorical ones use 0 for the same level and 1 otherwise.
- `<prefix>_mantel.csv` lists every Mantel r and p.

Mantel p values are not adjusted for the number of tests. With 3 targets × 19 variables, apply a correction (for example `p.adjust(..., "BH")` on the CSV) before reading single links too strictly.

## Changes in this revision

- Categorical variables were turned into arbitrary numeric codes. This put meaningless Pearson r values in the heatmap and fed code differences to the Mantel test. They are now left out of the heatmap and use a mismatch distance in the Mantel test.
- Line width was binned on r instead of |r|, so strong negative links were drawn as weak.
- The colour legend said "0.001 < p ≤ 0.01" for a class that also held p ≤ 0.001.
- Target node positions are computed for any number of targets, instead of being fixed for exactly Yield, PUE and SB.
- Input file and output prefix can be given on the command line, and missing columns give a clear error.
- The Mantel results are written to CSV.
- The README listed a helper script that is not in the repository. That entry is removed.

## Citation

```bibtex
@software{deng_mantel_cor_net,
  author = {Deng, Liangchao},
  title  = {mantel-cor\_net: correlation heatmap with Mantel-test links in R},
  year   = {2025},
  url    = {https://github.com/smiler488/mantel-cor_net}
}
```

## License

MIT, see [LICENSE](LICENSE). Liangchao Deng (邓良超).
