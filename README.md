# High-Dimensional Change-Point Detection

This repository contains the code and data-processing scripts used for the
numerical and empirical illustrations in the book chapter
*High-Dimensional Change-Point Detection* by Dan Zhuang and Long Feng.

## Repository contents

- `reproduce_illustrations.py`: Gaussian benchmark simulation and the CUSUM
  drift, power, and joint detection-localization figures.
- `reproduce_fred_md.R`: FRED-MD transformation, Double CUSUM analysis, result
  tables, and figure.
- `download_fred_md.R`: downloads the fixed December 2025 FRED-MD historical
  vintage from the official source and verifies its MD5 checksum.
- `reproduce_acgh.R`: ACGH global tests, multiple-change analyses, result
  tables, and figure.
- `plot_acgh_figure.R`: rebuilds the ACGH figure from saved analysis outputs.
- `results/`: saved summary results, metadata, random seeds, and software
  session information.
- `figures/`: figure outputs in publication and preview formats.
- `SOURCES_AND_LICENSES.md`: data provenance, method citations, third-party
  software licenses, and the boundary of this repository's MIT License.

The manuscript source is not part of this code repository.

## Software setup

The completed analyses used Python 3.12.3 and R 4.3.3 on Windows 11.

Install the Python dependencies in a virtual environment:

```sh
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
python -m pip install -r requirements.txt
```

Install the R dependencies, including the pinned development version of DMS:

```sh
Rscript install_r_dependencies.R
```

Exact package versions from the completed runs are recorded in
`results/simulation_metadata.json`, `results/fred_md_session_info.txt`, and
`results/acgh_session_info.txt`.

## Reproduce the results

Run all commands from the repository root.

### 1. Gaussian illustrations

```sh
python reproduce_illustrations.py
```

The script reuses completed simulation cells in `results/`. To rerun every
Monte Carlo cell, first move or remove the saved simulation result files.

### 2. FRED-MD example

Download the fixed official vintage and run the analysis:

```sh
Rscript download_fred_md.R
Rscript reproduce_fred_md.R
```

The first command downloads the official historical archive, extracts the
December 2025 monthly file, and saves it as `data/FRED-MD_2025m12.csv`.
See `FRED_MD_REPRODUCIBILITY.md` for the analysis specification and data-use
notes.

### 3. ACGH example

```sh
Rscript reproduce_acgh.R
```

The ACGH data are loaded from the `ecp` R package. See
`ACGH_REPRODUCIBILITY.md` for package versions, settings, and expected results.

The full bootstrap analyses can take substantially longer than rebuilding
figures from saved results. To regenerate only the ACGH figure, run:

```sh
Rscript plot_acgh_figure.R
```

## Reproducibility notes

- Random seeds and calibration settings are fixed in the scripts.
- The FRED-MD analysis uses a single bootstrap worker to preserve the random
  number stream.
- Saved results are included for comparison with a fresh run.
- Runtime can vary by platform; small Monte Carlo differences may occur when
  software versions or numerical libraries differ.

## Data and licensing

The source code is released under the MIT License. Third-party datasets are
not covered by that license.

FRED-MD data are obtained from the Federal Reserve Bank of St. Louis and may
contain series owned by third parties. The raw vintage and the transformed
panel are therefore not redistributed here. Users should review the FRED terms
and the terms attached to individual series before use. The ACGH data are
provided by the `ecp` R package and remain subject to the package and original
data-source terms.

See `SOURCES_AND_LICENSES.md` for the complete data and method provenance
record. No third-party source code is copied into this repository.

## Citation

Please cite the associated book chapter when its final bibliographic details
are available. GitHub can also generate a software citation from
`CITATION.cff`.

## Authors

- Dan Zhuang, School of Mathematics and Statistics, Fujian Normal University,
  Fuzhou 350117, China. Email: zhgdan@fjnu.edu.cn
- Long Feng, NITFID, School of Statistics and Data Science, LEBPS, LPMC, AAIS
  and KLMDASR, Nankai University, Tianjin, China. Email:
  flnankai@nankai.edu.cn
