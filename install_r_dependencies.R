#!/usr/bin/env Rscript

cran_versions <- c(
  ecp = "3.1.5",
  hdbinseg = "1.0.2",
  ggplot2 = "4.0.1",
  patchwork = "1.2.0",
  ragg = "1.3.0",
  digest = "0.6.37"
)

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}

for (package in names(cran_versions)) {
  target <- unname(cran_versions[[package]])
  installed <- requireNamespace(package, quietly = TRUE)
  current <- if (installed) as.character(utils::packageVersion(package)) else NA_character_
  if (!installed || current != target) {
    message("Installing ", package, " ", target)
    remotes::install_version(
      package,
      version = target,
      repos = "https://cloud.r-project.org",
      upgrade = "never"
    )
  }
}

dms_commit <- "79f4a9c07f48991943c5eb1cb9182e737041dfc2"
if (!requireNamespace("DMS", quietly = TRUE) ||
    as.character(utils::packageVersion("DMS")) != "0.0.0.9000") {
  message("Installing DMS from pinned GitHub commit ", dms_commit)
  remotes::install_github(
    paste0("ghwang-nk/DMS@", dms_commit),
    upgrade = "never"
  )
}

message("R dependencies are installed.")
