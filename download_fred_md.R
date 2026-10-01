#!/usr/bin/env Rscript

options(timeout = max(600, getOption("timeout")))

archive_url <- paste0(
  "https://www.stlouisfed.org/-/media/project/frbstl/stlouisfed/",
  "research/fred-md/historical-vintages-of-fred-md-2015-01-to-2025-12.zip"
)
target_file <- file.path("data", "FRED-MD_2025m12.csv")
expected_md5 <- "4626e1ab718adbfd623497e0aef5b191"

dir.create("data", showWarnings = FALSE, recursive = TRUE)

if (file.exists(target_file)) {
  observed_md5 <- unname(tools::md5sum(target_file))
  if (identical(observed_md5, expected_md5)) {
    message("Verified existing file: ", target_file)
    quit(status = 0L)
  }
  stop(
    "The existing data file has an unexpected MD5 checksum. Move or remove it ",
    "before downloading the fixed vintage."
  )
}

archive_file <- tempfile(fileext = ".zip")
extract_dir <- tempfile(pattern = "fred-md-")
dir.create(extract_dir)
on.exit(unlink(c(archive_file, extract_dir), recursive = TRUE, force = TRUE), add = TRUE)

message("Downloading the official FRED-MD historical-vintage archive...")
utils::download.file(archive_url, archive_file, mode = "wb", quiet = FALSE)

members <- utils::unzip(archive_file, list = TRUE)$Name
candidate <- members[grepl("(^|/)2025-12\\.csv$", members, ignore.case = TRUE)]
if (length(candidate) != 1L) {
  stop("Could not identify exactly one 2025-12.csv file in the archive.")
}

utils::unzip(archive_file, files = candidate, exdir = extract_dir)
source_file <- file.path(extract_dir, candidate)
if (!file.copy(source_file, target_file, overwrite = FALSE)) {
  stop("Failed to copy the extracted vintage to ", target_file)
}

observed_md5 <- unname(tools::md5sum(target_file))
if (!identical(observed_md5, expected_md5)) {
  unlink(target_file)
  stop(
    "Downloaded vintage failed the MD5 check. Expected ", expected_md5,
    "; observed ", observed_md5, "."
  )
}

message("Saved and verified: ", target_file)
