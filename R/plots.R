# Manuscript-style decomposition, pairwise, and driver plots.
#
# The three engine functions below (pace_pair_variance_pratt,
# plot_stacked_bar_focal_no_resp, make_pair_composite) are ported verbatim from
# the PACE analysis code (scripts/helpers/pace_pair_variance_pratt.R and
# scripts/core/06a-plots-decomp.R). Only the data plumbing changed: the exported
# wrappers read from a PACEFit rather than a bare `mv` bundle, so the figures
# reproduce the manuscript style exactly.

## ---------------------------------------------------------------------------
## Reconstruct the minimal `mv` bundle the ported functions expect.
## ---------------------------------------------------------------------------
.pace_mv_shim <- function(object) {
  gs <- colnames(object@fit$U)
  if (is.null(gs)) gs <- as.character(seq_len(ncol(object@fit$U)))
  list(fit = object@fit,
       gene_set = gs,
       mcsd_canonical = object@topDrivers,
       decomposition = list(
         gene_focal_single_frame = object@varianceDecomposition$perGene))
}

## SS-weighted per-focal block table (Goldstein pooling), from the single frame.
.pace_focal_blocks <- function(sf) {
  sf |>
    dplyr::mutate(ss_cell = .data[["Cell type %"]] / 100 * .data$denom,
                  ss_sp   = .data[["Spatial %"]]   / 100 * .data$denom,
                  ss_bl   = .data[["Spillover %"]] / 100 * .data$denom,
                  ss_rs   = .data[["Residual %"]]  / 100 * .data$denom) |>
    dplyr::group_by(.data$focal) |>
    dplyr::summarise(
      `Cell type`          = 100 * sum(.data$ss_cell) / sum(.data$denom),
      `Spatial cell state` = 100 * sum(.data$ss_sp)   / sum(.data$denom),
      `Spillover`          = 100 * sum(.data$ss_bl)   / sum(.data$denom),
      `Residual`           = 100 * sum(.data$ss_rs)   / sum(.data$denom),
      .groups = "drop")
}

## ---------------------------------------------------------------------------
## Pratt 1987 signed pair-level variance decomposition (verbatim engine port).
## ---------------------------------------------------------------------------
pace_pair_variance_pratt <- function(mv, cond_prefix = NULL, focals = NULL,
                                     cohort_label = "cohort", block_label = NULL) {
  fit <- mv$fit
  Z <- fit$re_meta$Z
  gn <- mv$gene_set
  colnames(fit$U) <- gn
  TYPES <- if (!is.null(fit$re_meta$blocks))
             fit$re_meta$blocks[[1]]$group_levels
           else fit$re_meta$group_levels
  if (is.null(focals)) focals <- TYPES
  if (is.null(block_label))
    block_label <- if (is.null(cond_prefix)) "Spatial" else "RxS"

  pair_rows <- list()
  focal_rows <- list()

  for (fc in focals) {
    fc_int <- paste0(fc, "::(Intercept)")
    if (!(fc_int %in% colnames(Z))) next
    cells <- which(as.numeric(Z[, fc_int]) != 0)
    if (length(cells) < 50) next
    term_names <- if (is.null(cond_prefix))
                    paste0(fc, "::", TYPES)
                  else paste0(fc, "::", cond_prefix, ":", TYPES)
    keep <- term_names %in% colnames(Z) & term_names %in% rownames(fit$U)
    if (!any(keep)) next
    tn <- term_names[keep]
    tt <- TYPES[keep]
    Z_fc <- as.matrix(Z[cells, tn, drop = FALSE])
    Sigma_K <- stats::cov(Z_fc)
    U_c <- as.matrix(fit$U[tn, , drop = FALSE]); U_c[!is.finite(U_c)] <- 0
    SU <- Sigma_K %*% U_c
    V_pair_gene <- U_c * SU
    V_pair_t <- as.numeric(rowSums(V_pair_gene))
    V_pratt_total <- sum(V_pair_t)
    V_diag_t <- as.numeric(rowSums(U_c^2) * diag(Sigma_K))
    V_diag_total <- sum(V_diag_t)

    for (i in seq_along(tt)) {
      pair_rows[[length(pair_rows) + 1]] <- data.frame(
        cohort = cohort_label, block = block_label,
        focal = fc, neighbour = tt[i],
        V_pair_pratt = V_pair_t[i],
        V_pair_diag  = V_diag_t[i],
        within_focal_share_pct = 100 * V_pair_t[i] / V_pratt_total,
        within_focal_diag_pct  = 100 * V_diag_t[i] / V_diag_total,
        sign = sign(V_pair_t[i]),
        stringsAsFactors = FALSE)
    }
    focal_rows[[length(focal_rows) + 1]] <- data.frame(
      cohort = cohort_label, block = block_label, focal = fc,
      n_cells = length(cells),
      V_block_pratt = V_pratt_total,
      V_block_diag  = V_diag_total,
      cross_cov_pct = 100 * (V_pratt_total - V_diag_total) / V_pratt_total,
      n_negative_pairs = sum(V_pair_t < 0),
      stringsAsFactors = FALSE)
  }
  list(pair_long = do.call(rbind, pair_rows),
       focal_summary = do.call(rbind, focal_rows))
}

## ---------------------------------------------------------------------------
## Per-focal stacked bar + zoom (verbatim engine port).
## ---------------------------------------------------------------------------
plot_stacked_bar_focal_no_resp <- function(
        blks_all,
        block_levels = c("Cell type", "Spatial cell state", "Spillover", "Residuals"),
        title = "") {
  keep_levels <- setdiff(block_levels, "Residuals")

  resid_by_focal <- blks_all %>%
    dplyr::filter(.data$block == "Residuals") %>%
    dplyr::mutate(focal = stringr::str_replace(.data$focal, "_", " ")) %>%
    dplyr::transmute(focal, resid_pct = .data$pct_total, explained = 1 - .data$pct_total)

  blks_plot <- blks_all %>%
    dplyr::filter(.data$block != "Residuals") %>%
    dplyr::mutate(focal = stringr::str_replace(.data$focal, "_", " "),
                  block = factor(.data$block, levels = keep_levels)) %>%
    tidyr::complete(focal, block, fill = list(pct_total = 0)) %>%
    dplyr::left_join(resid_by_focal, by = "focal")

  cols <- c("Cell type" = "#003049", "Spatial cell state" = "#780000",
            "Spillover" = "#E76F51")
  cols <- cols[keep_levels]

  focal_order <- blks_plot %>%
    dplyr::filter(block == "Spatial cell state") %>%
    dplyr::arrange(dplyr::desc(pct_total)) %>%
    dplyr::pull(focal)
  blks_plot$focal <- factor(blks_plot$focal, levels = focal_order)

  p_full <- ggplot(blks_plot, aes(x = focal, y = pct_total, fill = block)) +
    geom_col(width = 0.8, position = position_stack(reverse = TRUE)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = cols, breaks = keep_levels, limits = keep_levels, drop = FALSE) +
    labs(x = "Focal cell type", y = "Proportion of total variance",
         fill = "Block", title = title) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))

  zoom_data <- blks_plot %>%
    dplyr::filter(block %in% c("Spatial cell state", "Spillover"))

  p_zoom <- ggplot(zoom_data, aes(x = focal, y = pct_total, fill = block)) +
    geom_col(width = 0.7, position = position_dodge(width = 0.75)) +
    geom_text(aes(label = sprintf("%.2f%%", pct_total * 100)),
              position = position_dodge(width = 0.75), vjust = -0.4, size = 2.5) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                       expand = expansion(mult = c(0, 0.15))) +
    scale_fill_manual(values = cols, breaks = keep_levels, limits = keep_levels, drop = FALSE) +
    labs(x = "Focal cell type", y = NULL, fill = "Block",
         title = "Spatial cell state & Spillover") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 13),
          legend.position = "none")

  p_full <- p_full + theme(legend.position = "bottom")
  patchwork::wrap_plots(p_full, p_zoom, widths = c(1, 1))
}

## ---------------------------------------------------------------------------
## Per-pair MCSD driver bar + per-gene single-frame decomposition (verbatim).
## ---------------------------------------------------------------------------
make_pair_composite <- function(focus_focal, focus_neighbour, mv, sf = NULL, n_top = 5,
                                panels = c("both", "mcsd", "gene")) {
  panels      <- match.arg(panels)
  focus_key   <- paste(focus_focal, focus_neighbour, sep = "_")
  top_drivers <- utils::head(mv$mcsd_canonical[[focus_key]]$scores, n_top)

  p_mcsd <- ggplot2::ggplot(top_drivers, ggplot2::aes(MCSD, stats::reorder(gene, MCSD))) +
    ggplot2::geom_col(fill = "#F5A623", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", MCSD)), hjust = -0.15, size = 3) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.2))) +
    ggplot2::labs(x = "MCSD (spatial cell state)", y = NULL,
                  title = sprintf("%s <- %s", focus_focal, focus_neighbour)) +
    ggplot2::theme_classic(base_size = 11)
  if (panels == "mcsd") return(p_mcsd)

  if (is.null(sf)) stop("`sf` (single-frame decomposition) is required for the gene panel.")
  has_resp    <- "Responder spatial %" %in% names(sf)
  gene_blocks <- c("Cell type", "Spatial cell state",
                   if (has_resp) "Responder spatial state", "Spillover")
  block_cols  <- c("Cell type" = "#1F3B57", "Spatial cell state" = "#8B1A1A",
                   "Responder spatial state" = "#E07B39", "Spillover" = "#F6C9C9")

  gd <- dplyr::filter(sf, .data$focal == focus_focal, .data$gene %in% top_drivers$gene)
  gene_var <- data.frame(gene = gd$gene,
                         `Cell type`          = gd[["Cell type %"]],
                         `Spatial cell state` = gd[["Spatial %"]],
                         `Spillover`          = gd[["Spillover %"]],
                         check.names = FALSE)
  if (has_resp) gene_var[["Responder spatial state"]] <- gd[["Responder spatial %"]]
  gene_var <- gene_var |>
    tidyr::pivot_longer(-gene, names_to = "Block", values_to = "pct") |>
    dplyr::mutate(Block = factor(Block, levels = gene_blocks),
                  gene  = factor(gene, levels = top_drivers$gene))

  p_gene <- ggplot2::ggplot(gene_var, ggplot2::aes(gene, pct / 100, fill = Block)) +
    ggplot2::geom_col(width = 0.78, position = ggplot2::position_stack(reverse = TRUE)) +
    ggplot2::scale_fill_manual(values = block_cols[gene_blocks]) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(x = "Gene", y = "Proportion of total variance",
                  title = "Per-gene decomposition (single-frame)") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  if (panels == "gene") return(p_gene)

  patchwork::wrap_plots(p_mcsd, p_gene)
}

## ---------------------------------------------------------------------------
## Exported, PACEFit-facing wrappers
## ---------------------------------------------------------------------------

#' Per-focal variance decomposition plot
#'
#' Reproduces the manuscript per-focal decomposition figure: a full stacked bar
#' of the cell-type identity, spatial cell state, and contamination (spillover)
#' components, with a zoomed grouped bar of the two small spatial components.
#'
#' @param object A [PACEFit].
#' @param title Plot title.
#' @return A `patchwork` / `ggplot` object.
#' @examples
#' \dontrun{ plotDecomposition(fit) }
#' @export
plotDecomposition <- function(object, title = "Per-focal decomposition") {
  stopifnot(methods::is(object, "PACEFit"))
  sf <- object@varianceDecomposition$perGene
  fb <- .pace_focal_blocks(sf)
  blks_all <- fb |>
    dplyr::transmute(.data$focal,
                     `Cell type`          = .data$`Cell type` / 100,
                     `Spatial cell state` = .data$`Spatial cell state` / 100,
                     `Spillover`          = .data$Spillover / 100,
                     `Residuals`          = .data$Residual / 100) |>
    tidyr::pivot_longer(-.data$focal, names_to = "block", values_to = "pct_total")
  plot_stacked_bar_focal_no_resp(
    blks_all,
    block_levels = c("Cell type", "Spatial cell state", "Spillover", "Residuals"),
    title = title)
}

#' Pairwise spatial-percent heatmap
#'
#' Reproduces the manuscript pairwise heatmap: each cell is the percentage of the
#' focal type's total variance contributed by spatial interaction with that
#' neighbour, computed as the focal's spatial share split across neighbours by the
#' off-diagonal-renormalised Pratt attribution.
#'
#' @param object A [PACEFit].
#' @param title Plot title.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{ plotPairHeatmap(fit) }
#' @export
plotPairHeatmap <- function(object,
                            title = "Pairwise spatial % (variance, normalised Pratt)") {
  stopifnot(methods::is(object, "PACEFit"))
  mv <- .pace_mv_shim(object)
  sf <- object@varianceDecomposition$perGene
  fb <- .pace_focal_blocks(sf)
  TYPES  <- sort(unique(as.character(sf$focal)))
  foc_sp <- stats::setNames(fb$`Spatial cell state`, as.character(fb$focal))

  pratt_long <- pace_pair_variance_pratt(mv, cond_prefix = NULL,
                                         cohort_label = "cohort", block_label = "Spatial")$pair_long |>
    dplyr::filter(as.character(.data$focal) != as.character(.data$neighbour)) |>
    dplyr::group_by(.data$focal) |>
    dplyr::mutate(val = foc_sp[as.character(.data$focal)] *
                        .data$V_pair_pratt / sum(.data$V_pair_pratt)) |>
    dplyr::ungroup() |>
    dplyr::transmute(focal = as.character(.data$focal),
                     neighbour = as.character(.data$neighbour), .data$val)

  heat <- expand.grid(focal = TYPES, neighbour = TYPES, stringsAsFactors = FALSE) |>
    dplyr::filter(.data$focal != .data$neighbour) |>
    dplyr::left_join(pratt_long, by = c("focal", "neighbour")) |>
    dplyr::mutate(val = ifelse(is.na(.data$val), 0, .data$val),
                  focal = factor(.data$focal, levels = TYPES),
                  neighbour = factor(.data$neighbour, levels = rev(TYPES)))

  ggplot2::ggplot(heat, ggplot2::aes(.data$focal, .data$neighbour, fill = .data$val)) +
    ggplot2::geom_tile(colour = "grey92", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(abs(.data$val) < 0.005, "0.00",
                                                   sprintf("%.2f", .data$val))),
                       size = 2.4, colour = "white") +
    ggplot2::scale_fill_viridis_c(option = "inferno", na.value = "white",
                                  name = "spatial %\nof total") +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "Focal", y = "Neighbour", title = title) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.title = ggplot2::element_text(hjust = 0.5))
}

#' Per-pair driver composite
#'
#' Reproduces the manuscript per-pair driver figure: the top genes ranked by
#' driver score (MCSD) for a focal-neighbour pair, alongside their per-gene
#' single-frame variance decomposition.
#'
#' @param object A [PACEFit].
#' @param focal,neighbour The focal and neighbour cell types.
#' @param n_top Number of top genes to show.
#' @param panels One of `"both"`, `"mcsd"`, or `"gene"`.
#' @return A `patchwork` / `ggplot` object.
#' @examples
#' \dontrun{ plotDrivers(fit, "Macrophage", "Tumour") }
#' @export
plotDrivers <- function(object, focal, neighbour, n_top = 5,
                        panels = c("both", "mcsd", "gene")) {
  stopifnot(methods::is(object, "PACEFit"))
  mv <- .pace_mv_shim(object)
  sf <- object@varianceDecomposition$perGene
  make_pair_composite(focal, neighbour, mv, sf = sf, n_top = n_top,
                      panels = match.arg(panels))
}

## Canonical per-density-bin boxplot with the per-bin mean overlaid (verbatim
## style port of make_boxplot in build_bc_zoom_boxplots.R). For a zero-inflated
## gene the box median sits at 0, so the trend is carried by the mean (diamond +
## line); percent expressing per bin and the lm slope/p/n are annotated.
.pace_proximity_box <- function(d, gene_col, fill, ytitle, title) {
  dens_breaks <- c(-Inf, 0, 2, 5, 9, 14, 20, Inf)
  dens_labels <- c("0", "0-2", "2-5", "5-9", "9-14", "14-20", ">20")
  d$dens_bin <- cut(d$tumour_density, breaks = dens_breaks,
                    labels = dens_labels, right = TRUE)
  f  <- stats::lm(d[[gene_col]] ~ d$tumour_density)
  co <- summary(f)$coefficients
  slope <- co[2, 1]; pval <- co[2, 4]; n <- nrow(d)
  lev <- levels(d$dens_bin)
  pct <- tapply(d[[gene_col]], d$dens_bin, function(v) 100 * mean(v > 0))
  ypos <- min(d[[gene_col]]) - 0.06 * diff(range(d[[gene_col]]))
  lab_df  <- data.frame(dens_bin = factor(lev, levels = lev),
                        lab = sprintf("%.0f%%", pct[lev]), y = ypos)
  mean_df <- data.frame(dens_bin = factor(lev, levels = lev),
                        m = tapply(d[[gene_col]], d$dens_bin, mean)[lev])
  ann <- sprintf("slope = %+.3f\np = %s\nn = %d",
                 slope, format.pval(pval, digits = 2, eps = 1e-300), n)
  ggplot2::ggplot(d, ggplot2::aes(.data$dens_bin, .data[[gene_col]])) +
    ggplot2::geom_boxplot(fill = fill, colour = "grey25", alpha = 0.55,
                          outlier.size = 0.5, outlier.alpha = 0.3, linewidth = 0.4) +
    ggplot2::geom_line(data = mean_df, ggplot2::aes(.data$dens_bin, .data$m, group = 1),
                       inherit.aes = FALSE, colour = "black", linewidth = 0.7) +
    ggplot2::geom_point(data = mean_df, ggplot2::aes(.data$dens_bin, .data$m),
                        inherit.aes = FALSE, colour = "black", fill = "white",
                        shape = 23, size = 2.4, stroke = 0.8) +
    ggplot2::geom_text(data = lab_df, ggplot2::aes(.data$dens_bin, .data$y, label = .data$lab),
                       inherit.aes = FALSE, size = 2.5, colour = "grey30", vjust = 1) +
    ggplot2::annotate("label", x = 0.6, y = max(d[[gene_col]]), hjust = 0, vjust = 1,
                      label = ann, size = 2.9, fill = scales::alpha("white", 0.7)) +
    ggplot2::labs(
      x = "Tumour-neighbour density bin (low to high)  [% = focal cells expressing]",
      y = ytitle, title = paste0(title, "  (diamond = bin mean)")) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10),
                   axis.text.x = ggplot2::element_text(size = 6.5),
                   plot.margin = ggplot2::margin(6, 6, 14, 6))
}

#' Expression versus neighbour density (per-bin boxplot)
#'
#' Reproduces the manuscript proximity figure: for each focal cell, the
#' kernel-weighted density of a neighbour cell type is binned, and each gene's
#' log1p CP10k expression is shown as a boxplot per bin with the per-bin mean
#' overlaid (a diamond and connecting line), the percentage of focal cells
#' expressing the gene, and the fitted expression-vs-density slope.
#'
#' @param object A [PACEFit] (supplies the kernel bandwidths and cell-type
#'   column used at fit time).
#' @param spe The [SpatialExperiment::SpatialExperiment] that was fitted.
#' @param genes Character vector of genes to plot (one panel each).
#' @param focal,neighbour Focal and neighbour cell types.
#' @param assay_name Counts assay name (default `"counts"`).
#' @return A `patchwork` / `ggplot` object.
#' @examples
#' \dontrun{ plotProximity(fit, spe, c("MRC1", "APOC1"), "Macrophage", "Tumour") }
#' @export
plotProximity <- function(object, spe, genes, focal, neighbour,
                          assay_name = "counts") {
  stopifnot(methods::is(object, "PACEFit"))
  ct    <- as.character(SummarizedExperiment::colData(spe)[[object@params$celltype_col]])
  types <- sort(unique(ct))
  coords <- SpatialExperiment::spatialCoords(spe)
  h_bio <- object@params$h_bio; h_tech <- object@params$h_tech
  kern  <- pace_neighbour_kernel(coords, ct, types, h_bio = h_bio,
                                 h_tech = h_tech, eps = 3 * h_bio)
  dens  <- kern$K_bio[, which(types == neighbour)]

  counts <- SummarizedExperiment::assay(spe, assay_name)
  cp10k  <- log1p(Matrix::t(Matrix::t(counts) / pmax(Matrix::colSums(counts), 1) * 1e4))
  is_focal <- ct == focal
  fills <- c("#d62728", "#E76F51", "#1f77b4", "#2ca02c")

  plots <- lapply(seq_along(genes), function(i) {
    g <- genes[i]
    d <- data.frame(tumour_density = dens[is_focal],
                    val = as.numeric(cp10k[g, is_focal]))
    names(d)[2] <- g
    .pace_proximity_box(d, g, fills[(i - 1) %% length(fills) + 1],
                        sprintf("%s %s (log CP10k)", focal, g),
                        sprintf("%s in %s vs %s density", g, focal, neighbour))
  })
  patchwork::wrap_plots(plots, nrow = 1)
}
