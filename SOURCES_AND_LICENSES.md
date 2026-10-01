# Data Sources, Method Citations, and License Boundaries

This file records the provenance of every external dataset and named method
used in this repository. It also states which materials are and are not
covered by this repository's MIT License.

## Data sources

### FRED-MD macroeconomic panel

- **Provider:** Federal Reserve Bank of St. Louis, FRED-MD historical
  vintages.
- **Fixed input:** December 2025 monthly vintage. The analysis verifies the
  downloaded file against MD5 `4626e1ab718adbfd623497e0aef5b191`.
- **Official access:**
  https://www.stlouisfed.org/research/economists/mccracken/fred-databases
- **Dataset reference:** McCracken, M. W. and Ng, S. (2016), "FRED-MD: A
  Monthly Database for Macroeconomic Research," *Journal of Business &
  Economic Statistics*, 34(4), 574-589.
  https://doi.org/10.1080/07350015.2015.1086655
- **Rights notice:** FRED-MD combines series from multiple providers, and
  individual series can have different copyright and reuse conditions. The
  raw vintage and the transformed panel are not tracked in this repository.
  `download_fred_md.R` retrieves the vintage from the official source for the
  user. Users remain responsible for the FRED terms and series-specific
  restrictions: https://fred.stlouisfed.org/legal/terms/

### ACGH bladder-tumour panel

- **Immediate source:** the `ACGH` data object bundled with the R package
  `ecp` 3.1.5, loaded at run time by `data("ACGH", package = "ecp")`.
- **Package reference:** James, N. A. and Matteson, D. S. (2015), "ecp: An R
  Package for Nonparametric Multiple Change Point Analysis of Multivariate
  Data," *Journal of Statistical Software*, 62(7), 1-25.
  https://doi.org/10.18637/jss.v062.i07
- **Sources named by the `ecp::ACGH` help page:**
  - Bleakley, K. and Vert, J.-P. (2011), "The group fused Lasso for multiple
    change-point detection," arXiv:1106.4199.
    https://arxiv.org/abs/1106.4199
  - Stransky, N. et al. (2006), "Regional copy number-independent
    deregulation of transcription in cancer," *Nature Genetics*, 38,
    1386-1396. https://doi.org/10.1038/ng1923
- **Rights notice:** neither the complete ACGH matrix nor the patient-level
  figure-source CSV is tracked in this repository. They are regenerated
  locally from the installed `ecp` package. The `ecp` package reports license
  `GPL (>= 2)`; the original data and publications remain subject to their own
  terms.

## Method sources

### Gaussian benchmark illustrations

`reproduce_illustrations.py` is an independent implementation written for the
book chapter. It uses standard CUSUM contrasts with max and quadratic
aggregation and a Bonferroni combination. It is not copied from, and does not
claim to reimplement, DMS, INSPECT, or another named software method.

### FRED-MD transformations

The transformation codes in the FRED-MD file and their interpretation follow
McCracken and Ng (2016), cited above.

### DMS, max, sum, and adaptive tests

The ACGH script calls the separately installed `DMS` R package at pinned
commit `79f4a9c07f48991943c5eb1cb9182e737041dfc2`. The package is licensed
under GPL-3.0 and is not copied into this repository.

- **DMS, max, and sum procedures:** Wang, G. and Feng, L. (2023),
  "Computationally efficient and data-adaptive changepoint inference in high
  dimension," *Journal of the Royal Statistical Society Series B*, 85,
  936-958. https://doi.org/10.1093/jrsssb/qkad048
- **LZZL20 adaptive test:** Liu, B., Zhou, C., Zhang, X. and Liu, Y. (2020),
  "A Unified Data-Adaptive Framework for High Dimensional Change Point
  Detection," *Journal of the Royal Statistical Society Series B*, 82,
  933-963. https://doi.org/10.1111/rssb.12375
- **YC21 bootstrap max test:** Yu, M. and Chen, X. (2021), "Finite Sample
  Change Point Inference and Identification for High-Dimensional Mean
  Vectors," *Journal of the Royal Statistical Society Series B*, 83,
  247-270. https://doi.org/10.1111/rssb.12406

### Double CUSUM binary segmentation

Both empirical scripts call `hdbinseg::dcbs.alg`; the package is installed
separately and reports license `GPL (>= 3)`. The method source is:

- Cho, H. (2016), "Change-point detection in panel data via double CUSUM
  statistic," *Electronic Journal of Statistics*, 10(2), 2000-2038.
  https://doi.org/10.1214/16-EJS1155

## Software and license boundary

- The authors' scripts and documentation in this repository are released
  under the MIT License in `LICENSE`.
- Third-party R and Python packages are dependencies installed from their own
  repositories. They are not vendored here and remain under their own
  licenses.
- The MIT License does not relicense FRED-MD, ACGH, any third-party dataset,
  any third-party software, or any cited publication.
- Files under `figures/` and `results/` are supplied to verify the reported
  computations. Their inclusion does not grant rights in underlying
  third-party data.

This provenance record is intended to support responsible scholarly reuse; it
is not legal advice.
