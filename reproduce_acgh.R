#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 100)

# Data and method provenance: see SOURCES_AND_LICENSES.md.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "reproduce_acgh.R"
project_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
setwd(project_dir)

if (dir.exists(".Rlib")) {
  .libPaths(c(normalizePath(".Rlib", winslash = "/"), .libPaths()))
}

required_packages <- c("ecp", "DMS", "hdbinseg", "ggplot2", "patchwork", "ragg", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(
    "Missing packages: ", paste(missing_packages, collapse = ", "),
    ". See ACGH_REPRODUCIBILITY.md for installation instructions."
  )
}

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

data("ACGH", package = "ecp", envir = environment())
X <- ACGH$data
individual <- ACGH$individual
stopifnot(is.matrix(X), identical(dim(X), c(2215L, 43L)), length(individual) == ncol(X))
stopifnot(all(is.finite(X)))

n <- nrow(X)
p <- ncol(X)
boundary <- 40L
bootstrap_reps <- 500L
alpha_main <- 0.05
alpha_strict <- 0.01

panel <- data.frame(locus = seq_len(n), X, check.names = FALSE)
names(panel)[-1] <- sprintf("patient_%02d", individual)
write.csv(panel, "results/acgh_panel.csv", row.names = FALSE)

data_metadata <- data.frame(
  source_package = "ecp",
  source_package_version = as.character(utils::packageVersion("ecp")),
  source_package_url = "https://CRAN.R-project.org/package=ecp",
  source_package_license = packageDescription("ecp")[["License"]],
  object = "ACGH$data",
  package_citation = "James and Matteson (2015), doi:10.18637/jss.v062.i07",
  data_source_1 = "Bleakley and Vert (2011), arXiv:1106.4199",
  data_source_2 = "Stransky et al. (2006), doi:10.1038/ng1923",
  loci = n,
  profiles = p,
  patient_ids = paste(individual, collapse = ";"),
  serialized_matrix_sha256 = digest::digest(X, algo = "sha256", serialize = TRUE),
  stringsAsFactors = FALSE
)
write.csv(data_metadata, "results/acgh_data_metadata.csv", row.names = FALSE)

elapsed <- function(expr) {
  start <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, seconds = proc.time()[["elapsed"]] - start)
}

message("Running full-sequence tests ...")
dms0_time <- elapsed(DMS::DMS(X, gam = 0, type = "Fisher"))
dms05_time <- elapsed(DMS::DMS(X, gam = 0.5, lam = boundary, type = "Fisher"))
set.seed(20261001)
yc21_time <- elapsed(DMS::YC21(X, n.boundary = boundary, B = bootstrap_reps))

# Segment-specific seeds make a segment's bootstrap result independent of traversal order.
segment_seed <- function(start, end) {
  as.integer((20261002 + 1009 * start + 9176 * end) %% (.Machine$integer.max - 1L) + 1L)
}

lzzl_cache <- new.env(parent = emptyenv())
test_lzzl_segment <- function(start, end) {
  key <- paste(start, end, sep = ":")
  if (exists(key, envir = lzzl_cache, inherits = FALSE)) {
    return(get(key, envir = lzzl_cache, inherits = FALSE))
  }

  segment_length <- end - start + 1L
  if (segment_length < 2L * boundary) {
    result <- list(
      p_value = NA_real_, local_cp = NA_integer_, global_cp = NA_integer_,
      component_p_values = rep(NA_real_, 6L), seconds = 0
    )
  } else {
    set.seed(segment_seed(start, end))
    fit_time <- elapsed(DMS::LZZL20(
      X[start:end, , drop = FALSE],
      LZZL.s0 = p,
      LZZL.p = c(1:5, Inf),
      LZZL.B = bootstrap_reps,
      LZZL.boundary = boundary,
      est = TRUE
    ))
    fit <- fit_time$value
    local_cp <- as.integer(fit$th)
    result <- list(
      p_value = as.numeric(fit$pv),
      local_cp = local_cp,
      global_cp = start + local_cp - 1L,
      component_p_values = as.numeric(fit$pvs),
      seconds = fit_time$seconds
    )
  }
  assign(key, result, envir = lzzl_cache)
  result
}

run_lzzl_binary_segmentation <- function(alpha) {
  nodes <- list()
  changepoints <- integer()

  visit <- function(start, end, depth) {
    result <- test_lzzl_segment(start, end)
    rejected <- is.finite(result$p_value) && result$p_value <= alpha
    nodes[[length(nodes) + 1L]] <<- data.frame(
      alpha = alpha,
      start = start,
      end = end,
      length = end - start + 1L,
      depth = depth,
      p_value = result$p_value,
      changepoint = if (rejected) result$global_cp else NA_integer_,
      rejected = rejected,
      seconds = result$seconds
    )

    if (!rejected) return(invisible(NULL))
    cp <- result$global_cp
    if (!is.finite(cp) || cp < start + boundary - 1L || cp > end - boundary) {
      stop("Invalid LZZL20 split in interval [", start, ", ", end, "]: ", cp)
    }
    changepoints <<- c(changepoints, cp)
    visit(start, cp, depth + 1L)
    visit(cp + 1L, end, depth + 1L)
    invisible(NULL)
  }

  visit(1L, n, 0L)
  list(
    changepoints = sort(unique(changepoints)),
    nodes = do.call(rbind, nodes)
  )
}

message("Running adaptive LZZL20 binary segmentation ...")
lzzl_main <- run_lzzl_binary_segmentation(alpha_main)
lzzl_strict <- run_lzzl_binary_segmentation(alpha_strict)
root_lzzl <- test_lzzl_segment(1L, n)

message("Running Double CUSUM sensitivity analyses ...")
run_dcbs <- function(q, seed) {
  set.seed(seed)
  elapsed(hdbinseg::dcbs.alg(
    x = t(X),
    cp.type = 1,
    phi = -1,
    trim = boundary,
    height = 10,
    temporal = TRUE,
    B = 1000,
    q = q,
    do.parallel = 0
  ))
}

dcbs_main <- run_dcbs(alpha_main, 20261003)
dcbs_strict <- run_dcbs(alpha_strict, 20261004)
dcbs_main_cp <- sort(as.integer(dcbs_main$value$ecp))
dcbs_strict_cp <- sort(as.integer(dcbs_strict$value$ecp))

format_asymptotic_p <- function(x) {
  ifelse(x == 0, "<2.2e-16", formatC(x, digits = 4, format = "g"))
}
format_bootstrap_p <- function(x, B) {
  ifelse(x == 0, sprintf("<%.4f", 1 / (B + 1)), formatC(x, digits = 4, format = "f"))
}

global_tests <- data.frame(
  procedure = c(
    "Max(0)", "Sum", "DMS(0)", "Max(1/2)", "DMS(1/2)",
    "Yu--Chen Max(1/2)", "Liu--Zhou--Zhang--Liu adaptive"
  ),
  calibration = c(
    rep("asymptotic", 5),
    sprintf("multiplier bootstrap, B=%d", bootstrap_reps),
    sprintf("multiplier bootstrap, B=%d", bootstrap_reps)
  ),
  p_value = c(
    dms0_time$value$pv.Max,
    dms0_time$value$pv.Sum,
    dms0_time$value$pv.DMS,
    dms05_time$value$pv.Max,
    dms05_time$value$pv.DMS,
    yc21_time$value$pv,
    root_lzzl$p_value
  ),
  p_value_display = c(
    format_asymptotic_p(dms0_time$value$pv.Max),
    format_asymptotic_p(dms0_time$value$pv.Sum),
    format_asymptotic_p(dms0_time$value$pv.DMS),
    format_asymptotic_p(dms05_time$value$pv.Max),
    format_asymptotic_p(dms05_time$value$pv.DMS),
    format_bootstrap_p(yc21_time$value$pv, bootstrap_reps),
    format_bootstrap_p(root_lzzl$p_value, bootstrap_reps)
  ),
  seconds = c(
    dms0_time$seconds, dms0_time$seconds, dms0_time$seconds,
    dms05_time$seconds, dms05_time$seconds,
    yc21_time$seconds, root_lzzl$seconds
  ),
  stringsAsFactors = FALSE
)
write.csv(global_tests, "results/acgh_global_tests.csv", row.names = FALSE)

change_points <- rbind(
  data.frame(method = "LZZL20 adaptive binary segmentation", level = "alpha=0.05", locus = lzzl_main$changepoints),
  data.frame(method = "LZZL20 adaptive binary segmentation", level = "alpha=0.01", locus = lzzl_strict$changepoints),
  data.frame(method = "Double CUSUM combined", level = "q=0.05", locus = dcbs_main_cp),
  data.frame(method = "Double CUSUM combined", level = "q=0.01", locus = dcbs_strict_cp)
)
write.csv(change_points, "results/acgh_changepoints.csv", row.names = FALSE)
write.csv(lzzl_main$nodes, "results/acgh_lzzl_tree_alpha005.csv", row.names = FALSE)
write.csv(lzzl_strict$nodes, "results/acgh_lzzl_tree_alpha001.csv", row.names = FALSE)

nearest_distance <- function(query, reference) {
  if (!length(query)) return(integer())
  if (!length(reference)) return(rep(NA_integer_, length(query)))
  vapply(query, function(x) min(abs(reference - x)), integer(1))
}

agreement <- rbind(
  data.frame(
    method = "Double CUSUM combined", level = "q=0.05", locus = dcbs_main_cp,
    nearest_lzzl_distance = nearest_distance(dcbs_main_cp, lzzl_main$changepoints)
  ),
  data.frame(
    method = "Double CUSUM combined", level = "q=0.01", locus = dcbs_strict_cp,
    nearest_lzzl_distance = nearest_distance(dcbs_strict_cp, lzzl_main$changepoints)
  )
)
write.csv(agreement, "results/acgh_method_agreement.csv", row.names = FALSE)

summary_results <- data.frame(
  analysis = c("LZZL20 alpha=0.05", "LZZL20 alpha=0.01", "DCBS q=0.05", "DCBS q=0.01"),
  number_of_changepoints = c(
    length(lzzl_main$changepoints), length(lzzl_strict$changepoints),
    length(dcbs_main_cp), length(dcbs_strict_cp)
  ),
  changepoints = c(
    paste(lzzl_main$changepoints, collapse = ";"),
    paste(lzzl_strict$changepoints, collapse = ";"),
    paste(dcbs_main_cp, collapse = ";"),
    paste(dcbs_strict_cp, collapse = ";")
  ),
  seconds = c(
    sum(vapply(as.list(lzzl_cache), function(z) z$seconds, numeric(1))),
    0,
    dcbs_main$seconds,
    dcbs_strict$seconds
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_results, "results/acgh_analysis_summary.csv", row.names = FALSE)

message("Drawing the ACGH figure ...")
source("plot_acgh_figure.R", local = TRUE)
make_acgh_figure(project_dir = ".")

run_metadata <- data.frame(
  item = c(
    "run_time", "R_version", "platform", "ecp_version", "DMS_version",
    "hdbinseg_version", "ggplot2_version", "boundary_or_trim", "LZZL20_bootstrap_reps",
    "DCBS_bootstrap_reps", "LZZL20_s0", "LZZL20_norms", "DMS_gammas",
    "DMS_half_boundary", "DCBS_phi", "DCBS_height", "main_level", "strict_level"
  ),
  value = c(
    format(Sys.time(), tz = "UTC", usetz = TRUE),
    R.version.string,
    R.version$platform,
    as.character(utils::packageVersion("ecp")),
    as.character(utils::packageVersion("DMS")),
    as.character(utils::packageVersion("hdbinseg")),
    as.character(utils::packageVersion("ggplot2")),
    boundary,
    bootstrap_reps,
    1000,
    p,
    "1;2;3;4;5;Inf",
    "0;0.5",
    boundary,
    -1,
    10,
    alpha_main,
    alpha_strict
  ),
  stringsAsFactors = FALSE
)
write.csv(run_metadata, "results/acgh_run_metadata.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "results/acgh_session_info.txt")

message("Completed. Change-point counts:")
print(summary_results[, c("analysis", "number_of_changepoints")], row.names = FALSE)
message("Full-sequence p-values:")
print(global_tests[, c("procedure", "p_value_display")], row.names = FALSE)
