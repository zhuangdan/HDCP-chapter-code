# FRED-MD input data

The analysis uses the December 2025 monthly FRED-MD historical vintage.
The raw file is intentionally not redistributed in this repository because
individual FRED series may be subject to third-party terms.

From the repository root, run:

```sh
Rscript download_fred_md.R
```

The script downloads the official historical-vintage archive, extracts
`2025-12.csv`, saves it as `FRED-MD_2025m12.csv`, and verifies this MD5:

```text
4626e1ab718adbfd623497e0aef5b191
```

Official source:
https://www.stlouisfed.org/research/economists/mccracken/fred-databases

FRED terms:
https://fred.stlouisfed.org/legal/terms/
