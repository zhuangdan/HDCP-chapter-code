make_acgh_figure <- function(project_dir = ".") {
  old_working_directory <- setwd(project_dir)
  on.exit(setwd(old_working_directory), add = TRUE)

  if (dir.exists(".Rlib")) {
    .libPaths(c(normalizePath(".Rlib", winslash = "/"), .libPaths()))
  }

  required_packages <- c("ecp", "ggplot2", "patchwork", "ragg")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages)) {
    stop("Missing plotting packages: ", paste(missing_packages, collapse = ", "))
  }

  dir.create("results", showWarnings = FALSE, recursive = TRUE)
  dir.create("figures", showWarnings = FALSE, recursive = TRUE)

  data("ACGH", package = "ecp", envir = environment())
  X <- ACGH$data
  individual <- ACGH$individual
  n <- nrow(X)
  p <- ncol(X)
  stopifnot(identical(dim(X), c(2215L, 43L)))

  changepoints <- utils::read.csv("results/acgh_changepoints.csv", check.names = FALSE)
  agreement <- utils::read.csv("results/acgh_method_agreement.csv", check.names = FALSE)

  adaptive_cp <- sort(changepoints$locus[
    changepoints$method == "LZZL20 adaptive binary segmentation" &
      changepoints$level == "alpha=0.05"
  ])
  dcbs_main_cp <- sort(changepoints$locus[
    changepoints$method == "Double CUSUM combined" &
      changepoints$level == "q=0.05"
  ])
  dcbs_strict_cp <- sort(changepoints$locus[
    changepoints$method == "Double CUSUM combined" &
      changepoints$level == "q=0.01"
  ])
  stopifnot(length(adaptive_cp) == 42L, length(dcbs_main_cp) == 19L, length(dcbs_strict_cp) == 5L)

  # Four profiles are selected at equally spaced positions in file order, not by outcome.
  selected_columns <- unique(as.integer(round(seq(1, p, length.out = 4))))
  selected_ids <- individual[selected_columns]
  segment_starts <- c(1L, adaptive_cp + 1L)
  segment_ends <- c(adaptive_cp, n)

  fitted_profile <- function(values) {
    fitted <- numeric(length(values))
    for (segment in seq_along(segment_starts)) {
      index <- segment_starts[segment]:segment_ends[segment]
      fitted[index] <- stats::median(values[index])
    }
    fitted
  }

  profile_data <- do.call(rbind, lapply(seq_along(selected_columns), function(i) {
    column <- selected_columns[i]
    observed <- X[, column]
    data.frame(
      locus = seq_len(n),
      patient_id = selected_ids[i],
      patient = sprintf("Patient %d", selected_ids[i]),
      observed = observed,
      adaptive_segment_median = fitted_profile(observed),
      stringsAsFactors = FALSE
    )
  }))
  profile_data$patient <- factor(
    profile_data$patient,
    levels = sprintf("Patient %d", selected_ids)
  )
  write.csv(profile_data, "results/acgh_figure_profiles.csv", row.names = FALSE)

  agreement_main <- agreement[agreement$level == "q=0.05", , drop = FALSE]
  agreement_main$retained_at_1pct <- agreement_main$locus %in% dcbs_strict_cp
  agreement_main$status <- ifelse(
    agreement_main$retained_at_1pct,
    "Retained at 1%",
    "Detected at 5% only"
  )
  write.csv(agreement_main, "results/acgh_figure_agreement.csv", row.names = FALSE)

  strict_boundaries <- data.frame(locus = dcbs_strict_cp)
  clip_limit <- 1.75
  base_theme <- ggplot2::theme_classic(base_size = 8, base_family = "sans") +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.3, colour = "#222222"),
      axis.ticks = ggplot2::element_line(linewidth = 0.3, colour = "#222222"),
      axis.title = ggplot2::element_text(size = 8),
      axis.text = ggplot2::element_text(size = 7, colour = "#222222"),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 7),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 7, face = "bold", hjust = 0),
      panel.spacing.y = grid::unit(1.5, "mm"),
      plot.margin = ggplot2::margin(3, 5, 3, 3)
    )

  p_profiles <- ggplot2::ggplot(profile_data, ggplot2::aes(locus)) +
    ggplot2::geom_vline(
      data = strict_boundaries,
      ggplot2::aes(xintercept = locus, linetype = "Strict shared change point"),
      colour = "#D55E00", linewidth = 0.35, alpha = 0.8
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = observed, colour = "Observed log ratio"),
      linewidth = 0.18, alpha = 0.75
    ) +
    ggplot2::geom_step(
      ggplot2::aes(y = adaptive_segment_median, colour = "Adaptive segment median"),
      linewidth = 0.58, direction = "hv"
    ) +
    ggplot2::facet_wrap(~patient, ncol = 1, strip.position = "right") +
    ggplot2::scale_colour_manual(values = c(
      "Observed log ratio" = "#8E8E8E",
      "Adaptive segment median" = "#008B8B"
    )) +
    ggplot2::scale_linetype_manual(values = c("Strict shared change point" = "22")) +
    ggplot2::scale_x_continuous(
      limits = c(1, n), breaks = seq(0, 2000, 500), expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(breaks = c(-1, 0, 1)) +
    ggplot2::coord_cartesian(ylim = c(-clip_limit, clip_limit)) +
    ggplot2::labs(x = NULL, y = "Log ratio (clipped)", colour = NULL, linetype = NULL) +
    base_theme +
    ggplot2::theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      panel.grid.major.y = ggplot2::element_line(colour = "#E5E5E5", linewidth = 0.25),
      axis.line.x = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(order = 1, override.aes = list(linewidth = c(0.5, 0.8), alpha = 1)),
      linetype = ggplot2::guide_legend(order = 2, override.aes = list(colour = "#D55E00", linewidth = 0.5))
    )

  p_agreement <- ggplot2::ggplot(
    agreement_main,
    ggplot2::aes(locus, nearest_lzzl_distance)
  ) +
    ggplot2::annotate(
      "rect", xmin = 1, xmax = n, ymin = -0.7, ymax = 5,
      fill = "#E7F3F2", alpha = 0.65
    ) +
    ggplot2::geom_hline(yintercept = 0, colour = "#555555", linewidth = 0.3) +
    ggplot2::geom_hline(
      yintercept = 5, colour = "#008B8B", linetype = "22", linewidth = 0.35
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = locus, y = 0, yend = nearest_lzzl_distance),
      colour = "#B8B8B8", linewidth = 0.35
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = status), shape = 21, colour = "white",
      stroke = 0.35, size = 2.3
    ) +
    ggplot2::annotate(
      "text", x = 35, y = 19.5,
      label = "9 exact; 17 within 5 loci; all 19 within 20",
      hjust = 0, vjust = 1, size = 2.55, colour = "#333333"
    ) +
    ggplot2::scale_fill_manual(values = c(
      "Detected at 5% only" = "#666666",
      "Retained at 1%" = "#D55E00"
    )) +
    ggplot2::scale_x_continuous(
      limits = c(1, n), breaks = seq(0, 2000, 500), expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-0.7, 20), breaks = c(0, 5, 10, 15, 20),
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      x = "Ordered genomic locus",
      y = "Distance to nearest adaptive\nchange point (loci)",
      fill = NULL
    ) +
    base_theme +
    ggplot2::theme(
      legend.position = "top",
      panel.grid = ggplot2::element_blank()
    )

  figure <- p_profiles / p_agreement +
    patchwork::plot_layout(heights = c(3.2, 1.3)) +
    patchwork::plot_annotation(tag_levels = "a")

  ggplot2::ggsave(
    "figures/acgh_changepoints.pdf", figure,
    width = 183, height = 135, units = "mm", device = grDevices::cairo_pdf
  )
  ggplot2::ggsave(
    "figures/acgh_changepoints.svg", figure,
    width = 183, height = 135, units = "mm", device = grDevices::svg
  )
  ragg::agg_png(
    "figures/acgh_changepoints.png", width = 183, height = 135,
    units = "mm", res = 300, background = "white"
  )
  print(figure)
  grDevices::dev.off()
  ragg::agg_tiff(
    "figures/acgh_changepoints.tiff", width = 183, height = 135,
    units = "mm", res = 600, compression = "lzw", background = "white"
  )
  print(figure)
  grDevices::dev.off()

  invisible(list(
    figure = figure,
    selected_patient_ids = selected_ids,
    agreement = agreement_main
  ))
}

if (sys.nframe() == 0L) {
  make_acgh_figure(".")
}
