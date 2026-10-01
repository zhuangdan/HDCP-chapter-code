# FRED-MD reproducibility record

## Data

- Source: Federal Reserve Bank of St. Louis, FRED-MD historical vintages.
- Source page: https://www.stlouisfed.org/research/economists/mccracken/fred-databases
- Dataset reference: McCracken and Ng (2016), doi:10.1080/07350015.2015.1086655.
- Fixed vintage: December 2025 (`data/FRED-MD_2025m12.csv`). The raw file is
  not redistributed in this repository; obtain and verify it by running
  `Rscript download_fred_md.R`.
- Analysis window: January 1960 through August 2025.
- Missing-data rule: retain only series with finite officially transformed values throughout the window; no interpolation.
- Final dimensions: 788 months by 119 series.

## Run

From the repository root, run:

```text
Rscript download_fred_md.R
Rscript reproduce_fred_md.R
```

The script requires `hdbinseg`, `ggplot2`, `patchwork`, and `ragg`. The completed run used R 4.3.3 and `hdbinseg` 1.0.2. It uses fixed seeds and a single bootstrap worker so that the random-number stream is repeatable.

## Analysis

The script applies the official FRED-MD transformation codes, then runs Double CUSUM binary segmentation with `phi=0`, `phi=1/2`, and the combined rule (`phi=-1`). It separately tests means and diagonal second-order structure at Haar scale -1. The main bootstrap uses `B=1000` and `q=0.05`; a `q=0.01` sensitivity check is run for the detected `phi=1/2` second-order result.

The mean analysis detects no changes. For diagonal second-order structure, `phi=1/2` detects splits at 1979-08, 1982-08, 1985-03, 2020-01, and 2021-07. The combined rule detects 1979-08, 1982-08, and 1985-06; `phi=0` detects none. All five `phi=1/2` dates remain unchanged at `q=0.01`.

## Outputs

- `figures/fred_md_double_cusum.pdf`, `.png`, and `.svg`
- `results/fred_md_fit_summary.csv`
- `results/fred_md_changepoints.csv`
- `results/fred_md_changepoint_stability.csv`
- `results/fred_md_series_selection.csv`
- `results/fred_md_transformed_balanced_panel.csv`
- `results/fred_md_analysis_metadata.csv`
- `results/fred_md_session_info.txt`

Estimated dates are descriptive split locations, not confidence intervals or causal event dates. A date denotes the final month of the left segment.

The download helper retrieves the official historical-vintage archive from
the Federal Reserve Bank of St. Louis. FRED-MD series can have different
underlying data owners and terms. Users are responsible for reviewing the
FRED terms and any series-specific restrictions before reuse or
redistribution: https://fred.stlouisfed.org/legal/terms/

Full data and method citations, including the Double CUSUM source, are
collected in `SOURCES_AND_LICENSES.md`.
