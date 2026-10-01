#!/usr/bin/env Rscript

# Data and method provenance: see SOURCES_AND_LICENSES.md.

required <- c("hdbinseg", "ggplot2", "patchwork", "ragg")
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(hdbinseg)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE)
RNGkind("L'Ecuyer-CMRG")
base_seed <- 20260930L
bootstrap_B <- as.integer(Sys.getenv("FRED_MD_BOOTSTRAP_B", "1000"))
if (!is.finite(bootstrap_B) || bootstrap_B < 1L) {
  stop("FRED_MD_BOOTSTRAP_B must be a positive integer")
}
bootstrap_q <- as.numeric(Sys.getenv("FRED_MD_BOOTSTRAP_Q", "0.05"))
if (!is.finite(bootstrap_q) || bootstrap_q <= 0 || bootstrap_q >= 1) {
  stop("FRED_MD_BOOTSTRAP_Q must lie strictly between 0 and 1")
}

data_file <- file.path("data", "FRED-MD_2025m12.csv")
figure_dir <- "figures"
result_dir <- "results"
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(data_file)) {
  stop("Data file not found: ", data_file)
}

apply_tcode <- function(x, code) {
  n <- length(x)
  if (code == 1L) return(x)
  if (code == 2L) return(c(NA_real_, diff(x)))
  if (code == 3L) return(c(NA_real_, NA_real_, diff(x, differences = 2L)))
  if (code %in% 4:6 && any(x[is.finite(x)] <= 0)) {
    stop("Non-positive value encountered for logarithmic transformation")
  }
  if (code == 4L) return(log(x))
  if (code == 5L) return(c(NA_real_, diff(log(x))))
  if (code == 6L) return(c(NA_real_, NA_real_, diff(log(x), differences = 2L)))
  if (code == 7L) {
    growth <- c(NA_real_, x[-1L] / x[-n] - 1)
    return(c(NA_real_, diff(growth)))
  }
  stop("Unknown FRED-MD transformation code: ", code)
}

raw_file <- read.csv(
  data_file,
  check.names = FALSE,
  na.strings = c("", "NA", "#N/A"),
  stringsAsFactors = FALSE
)

series_names <- names(raw_file)[-1L]
tcodes <- as.integer(unlist(raw_file[1L, -1L], use.names = FALSE))
if (length(tcodes) != length(series_names) || any(!tcodes %in% 1:7)) {
  stop("Invalid transformation-code row in ", data_file)
}

dates_all <- as.Date(raw_file[-1L, 1L], format = "%m/%d/%Y")
levels_all <- as.data.frame(
  lapply(raw_file[-1L, -1L, drop = FALSE], as.numeric),
  check.names = FALSE
)
transformed_all <- vapply(
  seq_along(series_names),
  function(j) apply_tcode(levels_all[[j]], tcodes[j]),
  numeric(length(dates_all))
)
colnames(transformed_all) <- series_names

analysis_start <- as.Date("1960-01-01")
analysis_end <- as.Date("2025-08-01")
window <- dates_all >= analysis_start & dates_all <= analysis_end
dates <- dates_all[window]
transformed_window <- transformed_all[window, , drop = FALSE]

finite_series <- apply(transformed_window, 2L, function(x) all(is.finite(x)))
nonconstant_series <- apply(
  transformed_window[, finite_series, drop = FALSE],
  2L,
  function(x) is.finite(sd(x)) && sd(x) > 0
)
keep <- finite_series
keep[finite_series] <- nonconstant_series

panel <- transformed_window[, keep, drop = FALSE]
if (ncol(panel) < 2L || nrow(panel) < 3L) {
  stop("The balanced transformed panel is too small")
}

# hdbinseg expects component series in rows and time in columns.
x <- t(panel)
method_spec <- data.frame(
  method = c("phi=0", "phi=1/2", "combined"),
  phi = c(0, 0.5, -1),
  seed = base_seed + 0:2,
  stringsAsFactors = FALSE
)

run_dcbs <- function(feature, q, seed_offset) {
  cp_type <- if (feature == "mean") 1L else 2L
  fits <- vector("list", nrow(method_spec))
  for (i in seq_len(nrow(method_spec))) {
    set.seed(method_spec$seed[i] + seed_offset)
    message(
      "Running Double CUSUM: feature=", feature,
      ", method=", method_spec$method[i],
      ", q=", q, " (B=", bootstrap_B, ")"
    )
    arguments <- list(
      x = x,
      cp.type = cp_type,
      phi = method_spec$phi[i],
      temporal = TRUE,
      B = bootstrap_B,
      q = q,
      do.parallel = 0
    )
    if (cp_type == 2L) {
      arguments$scales <- -1L
      arguments$diag <- TRUE
    }
    fits[[i]] <- do.call(hdbinseg::dcbs.alg, arguments)
  }
  names(fits) <- method_spec$method
  fits
}

summarize_fits <- function(fits, feature, q, analysis) {
  data.frame(
    analysis = analysis,
    feature = feature,
    method = method_spec$method,
    phi = method_spec$phi,
    bootstrap_q = q,
    bootstrap_threshold = vapply(fits, function(fit) as.numeric(fit$thr), numeric(1)),
    detected_changepoints = vapply(fits, function(fit) length(fit$ecp), integer(1)),
    stringsAsFactors = FALSE
  )
}

extract_changepoints <- function(fits, feature, q, analysis) {
  rows <- lapply(seq_along(fits), function(i) {
    ecp <- as.integer(fits[[i]]$ecp)
    if (length(ecp) == 0L) return(NULL)
    data.frame(
      analysis = analysis,
      feature = feature,
      method = method_spec$method[i],
      phi = method_spec$phi[i],
      bootstrap_q = q,
      split_index = ecp,
      left_segment_end = dates[ecp],
      right_segment_start = dates[pmin(ecp + 1L, length(dates))],
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    return(data.frame(
      analysis = character(), feature = character(), method = character(),
      phi = numeric(), bootstrap_q = numeric(), split_index = integer(),
      left_segment_end = as.Date(character()),
      right_segment_start = as.Date(character()),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

mean_fits <- run_dcbs("mean", bootstrap_q, 0L)
second_fits <- run_dcbs("diagonal_second_order", bootstrap_q, 100L)

# A stricter threshold checks whether the dense second-order dates depend on q.
strict_spec <- method_spec[method_spec$method == "phi=1/2", , drop = FALSE]
set.seed(base_seed + 200L)
message("Running strict-threshold sensitivity: feature=diagonal_second_order, method=phi=1/2, q=0.01")
strict_second_fit <- hdbinseg::dcbs.alg(
  x,
  cp.type = 2,
  phi = strict_spec$phi,
  scales = -1L,
  diag = TRUE,
  B = bootstrap_B,
  q = 0.01,
  do.parallel = 0
)

fit_summary <- rbind(
  summarize_fits(mean_fits, "mean", bootstrap_q, "main"),
  summarize_fits(second_fits, "diagonal_second_order", bootstrap_q, "main"),
  data.frame(
    analysis = "strict_threshold_sensitivity",
    feature = "diagonal_second_order",
    method = "phi=1/2",
    phi = 0.5,
    bootstrap_q = 0.01,
    bootstrap_threshold = as.numeric(strict_second_fit$thr),
    detected_changepoints = length(strict_second_fit$ecp),
    stringsAsFactors = FALSE
  )
)
write.csv(fit_summary, file.path(result_dir, "fred_md_fit_summary.csv"), row.names = FALSE)

changepoints <- rbind(
  extract_changepoints(mean_fits, "mean", bootstrap_q, "main"),
  extract_changepoints(second_fits, "diagonal_second_order", bootstrap_q, "main"),
  data.frame(
    analysis = "strict_threshold_sensitivity",
    feature = "diagonal_second_order",
    method = "phi=1/2",
    phi = 0.5,
    bootstrap_q = 0.01,
    split_index = as.integer(strict_second_fit$ecp),
    left_segment_end = dates[as.integer(strict_second_fit$ecp)],
    right_segment_start = dates[pmin(as.integer(strict_second_fit$ecp) + 1L, length(dates))],
    stringsAsFactors = FALSE
  )
)
write.csv(changepoints, file.path(result_dir, "fred_md_changepoints.csv"), row.names = FALSE)

main_dense <- changepoints[
  changepoints$analysis == "main" & changepoints$feature == "diagonal_second_order" &
    changepoints$method == "phi=1/2",
  , drop = FALSE
]
strict_dense <- changepoints[
  changepoints$analysis == "strict_threshold_sensitivity",
  , drop = FALSE
]
stability <- merge(
  main_dense[, c("split_index", "left_segment_end"), drop = FALSE],
  strict_dense[, c("split_index", "left_segment_end"), drop = FALSE],
  by = "split_index", all = TRUE, suffixes = c("_q005", "_q001")
)
stability$detected_at_both_thresholds <-
  !is.na(stability$left_segment_end_q005) & !is.na(stability$left_segment_end_q001)
write.csv(
  stability,
  file.path(result_dir, "fred_md_changepoint_stability.csv"),
  row.names = FALSE
)

selection <- data.frame(
  series = series_names,
  transformation_code = tcodes,
  retained = keep,
  stringsAsFactors = FALSE
)
write.csv(selection, file.path(result_dir, "fred_md_series_selection.csv"), row.names = FALSE)

processed_panel <- data.frame(date = dates, panel, check.names = FALSE)
write.csv(
  processed_panel,
  file.path(result_dir, "fred_md_transformed_balanced_panel.csv"),
  row.names = FALSE
)

metadata <- data.frame(
  field = c(
    "data_source", "source_page", "source_terms", "dataset_citation",
    "vintage", "source_file", "source_md5", "analysis_start",
    "analysis_end", "n_months", "p_series", "transformation", "missing_data_rule",
    "method", "mean_feature", "second_order_feature", "hdbinseg_version",
    "R_version", "bootstrap_B", "bootstrap_q", "sensitivity_q",
    "temporal_scaling", "parallel_workers", "base_seed", "default_trim",
    "default_tree_height"
  ),
  value = c(
    "Federal Reserve Bank of St. Louis FRED-MD historical vintages",
    "https://www.stlouisfed.org/research/economists/mccracken/fred-databases",
    "https://fred.stlouisfed.org/legal/terms/",
    "McCracken and Ng (2016), doi:10.1080/07350015.2015.1086655",
    "2025-12", basename(data_file), unname(tools::md5sum(data_file)),
    as.character(analysis_start), as.character(analysis_end), nrow(panel), ncol(panel),
    "Official FRED-MD transformation codes 1-7",
    "Retain only series finite over the prespecified window; no interpolation",
    "Double CUSUM binary segmentation",
    "Changes in the means of officially transformed series",
    "Changes in diagonal second-order structure at Haar scale -1",
    as.character(utils::packageVersion("hdbinseg")), R.version.string,
    bootstrap_B, bootstrap_q, 0.01, TRUE, 0, base_seed,
    round(log(nrow(panel))^1.5), floor(log(nrow(panel), 2) / 2)
  ),
  stringsAsFactors = FALSE
)
write.csv(metadata, file.path(result_dir, "fred_md_analysis_metadata.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(result_dir, "fred_md_session_info.txt"))

z_panel <- abs(scale(panel))
z_panel[z_panel > 3] <- 3
p_series <- ncol(z_panel)
heatmap_data <- data.frame(
  time_index = rep(seq_along(dates), times = p_series),
  series_index = rep(seq_len(p_series), each = length(dates)),
  magnitude = as.vector(z_panel),
  stringsAsFactors = FALSE
)

y_breaks <- unique(round(seq(1, p_series, length.out = min(9L, p_series))))
primary_indices <- main_dense$split_index
time_breaks <- unique(round(seq(1, length(dates), length.out = 8L)))

palette <- c(
  blue = "#2B6F9E",
  pale = "#F7F7F4",
  orange = "#D47A2C",
  dark = "#252525",
  grey = "#777777"
)

base_theme <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = palette["dark"]),
    axis.ticks = element_line(linewidth = 0.35, colour = palette["dark"]),
    axis.text = element_text(colour = palette["dark"]),
    axis.title = element_text(colour = palette["dark"]),
    plot.tag = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 7.5),
    legend.text = element_text(size = 7),
    plot.margin = margin(3, 5, 3, 3)
  )

p_heat <- ggplot(heatmap_data, aes(x = time_index, y = series_index, fill = magnitude)) +
  geom_raster(interpolate = FALSE) +
  geom_vline(
    xintercept = primary_indices,
    colour = palette["dark"], linewidth = 0.35, alpha = 0.8
  ) +
  scale_fill_gradient(
    low = palette["pale"], high = palette["blue"],
    limits = c(0, 3), name = "Absolute standardized\ntransformed value"
  ) +
  scale_y_continuous(
    breaks = y_breaks,
    labels = colnames(panel)[y_breaks],
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_continuous(
    breaks = time_breaks, labels = format(dates[time_breaks], "%Y"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.line.x = element_blank(), legend.position = "right"
  )

method_levels <- c("phi=0", "phi=1/2", "combined")
timeline_data <- changepoints[
  changepoints$analysis == "main" & changepoints$feature == "diagonal_second_order",
  , drop = FALSE
]
timeline_data$method <- factor(timeline_data$method, levels = rev(method_levels))
p_timeline <- ggplot(timeline_data, aes(x = split_index, y = method)) +
  geom_hline(
    yintercept = seq_along(method_levels),
    colour = "#DDDDDD", linewidth = 0.3
  ) +
  geom_point(
    aes(shape = method, colour = method),
    size = 2.1, stroke = 0.45
  ) +
  scale_shape_manual(
    values = c("combined" = 15, "phi=1/2" = 16, "phi=0" = 17),
    guide = "none"
  ) +
  scale_colour_manual(
    values = c(
      "combined" = unname(palette["dark"]),
      "phi=1/2" = unname(palette["orange"]),
      "phi=0" = unname(palette["blue"])
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(1, length(dates)), breaks = time_breaks,
    labels = format(dates[time_breaks], "%Y"),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(limits = rev(method_levels), drop = FALSE) +
  labs(x = "Month", y = NULL) +
  base_theme +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

figure <- p_heat / p_timeline +
  plot_layout(heights = c(4.2, 1)) +
  plot_annotation(tag_levels = "a")

figure_base <- file.path(figure_dir, "fred_md_double_cusum")
ggsave(
  paste0(figure_base, ".pdf"), figure,
  device = grDevices::cairo_pdf, width = 7.2, height = 5.9, units = "in",
  family = "Arial"
)
ggsave(
  paste0(figure_base, ".png"), figure,
  device = ragg::agg_png, width = 7.2, height = 5.9, units = "in", dpi = 300
)
ggsave(
  paste0(figure_base, ".svg"), figure,
  device = grDevices::svg, width = 7.2, height = 5.9, units = "in",
  family = "Arial"
)

cat("Analysis complete\n")
cat("Window:", as.character(min(dates)), "to", as.character(max(dates)), "\n")
cat("Panel dimensions: n =", nrow(panel), ", p =", ncol(panel), "\n")
for (feature in c("mean", "diagonal_second_order")) {
  for (i in seq_len(nrow(method_spec))) {
    cp_dates <- changepoints$left_segment_end[
      changepoints$analysis == "main" & changepoints$feature == feature &
        changepoints$method == method_spec$method[i]
    ]
    cat(feature, ", ", method_spec$method[i], ": ", paste(cp_dates, collapse = ", "), "\n", sep = "")
  }
}
