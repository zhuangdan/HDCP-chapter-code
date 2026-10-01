# Reproducing the ACGH example

The chapter analysis uses the `ACGH` data bundled with `ecp` 3.1.5. The object contains 2,215 ordered genomic loci for 43 bladder-tumour profiles. No observations are deleted, imputed, standardized, or reordered. The `ecp::ACGH` help page attributes the data to Bleakley and Vert (2011, arXiv:1106.4199) and Stransky et al. (2006, doi:10.1038/ng1923); the package itself is described by James and Matteson (2015, doi:10.18637/jss.v062.i07). Full data and method citations are collected in `SOURCES_AND_LICENSES.md`.

## Software

The analysis was run with R 4.3.3 and the following packages:

- `ecp` 3.1.5
- `DMS` 0.0.0.9000, GitHub commit `79f4a9c07f48991943c5eb1cb9182e737041dfc2`
- `hdbinseg` 1.0.2
- `ggplot2`, `patchwork`, `ragg`, and `digest`

Install the CRAN packages and the pinned development version of `DMS` with:

```r
install.packages(c("ecp", "hdbinseg", "ggplot2", "patchwork", "ragg", "digest", "remotes"))
remotes::install_github("ghwang-nk/DMS@79f4a9c07f48991943c5eb1cb9182e737041dfc2")
```

## Run

From this directory, run:

```sh
Rscript reproduce_acgh.R
```

The script verifies the expected dimensions and finite values, writes an auditable local CSV copy of the analysed matrix, runs all tests and segmentations, and recreates the figure in PDF, SVG, PNG, and TIFF formats. The matrix and patient-level figure-source CSV are ignored by Git and are not redistributed. Seeds are fixed. LZZL20 uses segment-specific seeds, so a segment's multiplier-bootstrap result does not depend on recursion order.

## Prespecified analyses

The full-sequence table contains the max-, sum-, and double-max-sum tests in `DMS`, the bootstrap Max(1/2) test of Yu and Chen, and the adaptive test of Liu, Zhou, Zhang, and Liu. For the latter, the implementation follows the published ACGH analysis: `s0 = 43`, norms `1, 2, 3, 4, 5, Inf`, 500 multiplier-bootstrap samples, and a 40-locus boundary.

Multiple changes are estimated by recursively applying that adaptive test at levels 0.05 and 0.01. A Double CUSUM comparison uses the combined weighting rule (`phi = -1`), mean-change target, temporal scale estimates, 1,000 bootstrap samples, a 40-locus trim, and the same two tail probabilities. The Double CUSUM tree height is set to 10 so that the software default does not impose an artificial cap below the number of published ACGH change points.

The methods have different statistics and calibration theories. Agreement is therefore descriptive, not a measure of truth. The data object has ordered locus indices but no chromosome or gene annotation, so estimated indices are not assigned biological labels.

With the pinned software and seeds, all full-sequence tests reject. Adaptive binary segmentation returns 42 change points at both levels, compared with 43 in the original published analysis whose random seed was not reported. Combined Double CUSUM returns 19 and 5 change points at the 0.05 and 0.01 thresholds, respectively. All five strict Double CUSUM locations exactly match adaptive locations.

## Outputs

- `results/acgh_data_metadata.csv`: data provenance, dimensions, patient IDs, and matrix hash
- `results/acgh_panel.csv`: locally generated analysed 2,215 by 43 matrix (not tracked)
- `results/acgh_global_tests.csv`: full-sequence p-values and runtimes
- `results/acgh_changepoints.csv`: all reported change-point indices
- `results/acgh_lzzl_tree_alpha005.csv` and `results/acgh_lzzl_tree_alpha001.csv`: tested recursive intervals
- `results/acgh_method_agreement.csv`: nearest-method distances for Double CUSUM change points
- `results/acgh_figure_profiles.csv`: locally generated patient-level figure data (not tracked)
- `results/acgh_figure_agreement.csv`: derived method-agreement data for the figure
- `results/acgh_analysis_summary.csv`: counts and complete change-point lists
- `results/acgh_run_metadata.csv` and `results/acgh_session_info.txt`: settings and software environment
- `figures/acgh_changepoints.*`: publication and preview versions of the figure

The figure can be regenerated from existing result tables without rerunning the bootstrap calculations by running `Rscript plot_acgh_figure.R`.
