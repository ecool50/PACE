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
  has_resp <- "Responder spatial %" %in% names(sf)
  fb <- sf |>
    dplyr::mutate(
      ss_cell = .data[["Cell type %"]] / 100 * .data$denom,
      ss_sp   = .data[["Spatial %"]]   / 100 * .data$denom,
      ss_bl   = .data[["Spillover %"]] / 100 * .data$denom,
      ss_rs   = .data[["Residual %"]]  / 100 * .data$denom,
      ss_rsp  = if (has_resp) .data[["Responder spatial %"]] / 100 * .data$denom
                else 0) |>
    dplyr::group_by(.data$focal) |>
    dplyr::summarise(
      `Cell type`               = 100 * sum(.data$ss_cell) / sum(.data$denom),
      `Spatial cell state`      = 100 * sum(.data$ss_sp)   / sum(.data$denom),
      `Responder spatial state` = 100 * sum(.data$ss_rsp)  / sum(.data$denom),
      `Spillover`               = 100 * sum(.data$ss_bl)   / sum(.data$denom),
      `Residual`                = 100 * sum(.data$ss_rs)   / sum(.data$denom),
      .groups = "drop")
  if (!has_resp) fb[["Responder spatial state"]] <- NULL
  fb
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
## Per-focal stacked bar + zoom WITH the responder spatial block (verbatim port
## of plot_stacked_bar_focal, for condition cohorts).
## ---------------------------------------------------------------------------
plot_stacked_bar_focal <- function(
        blks_all,
        block_levels = c("Cell type", "Spatial cell state", "Responder spatial state",
                         "Responder status", "Spillover", "Residuals"),
        title = "",
        zoom_blocks = c("Spatial cell state", "Responder spatial state", "Spillover")) {
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
            "Responder spatial state" = "#c1121f", "Responder status" = "#669bbc",
            "Spillover" = "#E5A100")
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

  zoom_blocks <- intersect(zoom_blocks, keep_levels)
  zoom_data <- blks_plot %>% dplyr::filter(block %in% zoom_blocks)

  p_zoom <- ggplot(zoom_data, aes(x = focal, y = pct_total, fill = block)) +
    geom_col(width = 0.7, position = position_dodge(width = 0.75)) +
    geom_text(aes(label = sprintf("%.2f%%", pct_total * 100)),
              position = position_dodge(width = 0.75), vjust = -0.4, size = 2.5) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                       expand = expansion(mult = c(0, 0.15))) +
    scale_fill_manual(values = cols, breaks = keep_levels, limits = keep_levels, drop = FALSE) +
    labs(x = "Focal cell type", y = NULL, fill = "Block",
         title = paste(zoom_blocks, collapse = " & ")) +
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
  has_resp <- "Responder spatial state" %in% names(fb)

  if (has_resp) {
    blks_all <- fb |>
      dplyr::transmute(.data$focal,
                       `Cell type`               = .data$`Cell type` / 100,
                       `Spatial cell state`      = .data$`Spatial cell state` / 100,
                       `Responder spatial state` = .data$`Responder spatial state` / 100,
                       `Spillover`               = .data$Spillover / 100,
                       `Residuals`               = .data$Residual / 100) |>
      tidyr::pivot_longer(-.data$focal, names_to = "block", values_to = "pct_total")
    plot_stacked_bar_focal(
      blks_all,
      block_levels = c("Cell type", "Spatial cell state", "Responder spatial state",
                       "Spillover", "Residuals"),
      title = title)
  } else {
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
}

#' Pairwise spatial-percent heatmap
#'
#' Reproduces the manuscript pairwise heatmap: each cell is the percentage of the
#' focal type's total variance contributed by spatial interaction with that
#' neighbour, computed as the focal's spatial share split across neighbours by the
#' off-diagonal-renormalised Pratt attribution.
#'
#' @param object A [PACEFit].
#' @param block Which spatial block to map: `"spatial"` (default, the baseline
#'   spatial cell state) or `"responder"` (the responder-by-proximity block, for
#'   condition cohorts).
#' @param title Plot title; a sensible default is chosen per `block`.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{ plotPairHeatmap(fit); plotPairHeatmap(fit, block = "responder") }
#' @export
plotPairHeatmap <- function(object, block = c("spatial", "responder"), title = NULL) {
  stopifnot(methods::is(object, "PACEFit"))
  block <- match.arg(block)
  mv <- .pace_mv_shim(object)
  sf <- object@varianceDecomposition$perGene
  fb <- .pace_focal_blocks(sf)
  TYPES <- sort(unique(as.character(sf$focal)))

  if (block == "responder") {
    if (!"Responder spatial state" %in% names(fb))
      stop("this fit has no responder block; use block = 'spatial'.", call. = FALSE)
    cond_prefix <- object@params$resp_term
    foc         <- stats::setNames(fb$`Responder spatial state`, as.character(fb$focal))
    block_label <- "RxS"
    fill_name   <- "responder %\nof total"
    if (is.null(title)) title <- "Pairwise responder spatial % (normalised Pratt)"
  } else {
    cond_prefix <- NULL
    foc         <- stats::setNames(fb$`Spatial cell state`, as.character(fb$focal))
    block_label <- "Spatial"
    fill_name   <- "spatial %\nof total"
    if (is.null(title)) title <- "Pairwise spatial % (variance, normalised Pratt)"
  }

  pratt_long <- pace_pair_variance_pratt(mv, cond_prefix = cond_prefix,
                                         cohort_label = "cohort",
                                         block_label = block_label)$pair_long |>
    dplyr::filter(as.character(.data$focal) != as.character(.data$neighbour)) |>
    dplyr::group_by(.data$focal) |>
    dplyr::mutate(val = foc[as.character(.data$focal)] *
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
                                  name = fill_name) +
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

#' Expression versus number of neighbours (per-bin boxplot)
#'
#' Reproduces the manuscript proximity figure (verbatim style port of
#' `density_bin_boxplot` in the analysis helpers): for each focal cell the number
#' of `neighbour`-type cells within `radius` micrometres is counted and binned,
#' and each gene's raw counts are drawn as a boxplot per bin with the per-bin mean
#' overlaid as a point. Because these genes are zero-inflated the box median sits
#' at zero in most bins, so the mean point carries the trend.
#'
#' @param object A [PACEFit] (supplies the cell-type column used at fit time).
#' @param spe The [SpatialExperiment::SpatialExperiment] that was fitted.
#' @param genes Character vector of genes (shown as facets).
#' @param focal,neighbour Focal and neighbour cell types.
#' @param radius Neighbour search radius in micrometres (default 30).
#' @param breaks Neighbour-count bin breaks (default `c(0, 2, 4, 6, 8, Inf)`).
#' @param box_colour Box and mean-point colour.
#' @param assay_name Counts assay name (default `"counts"`).
#' @return A `ggplot` object.
#' @examples
#' \dontrun{ plotProximity(fit, spe, c("MRC1", "APOC1"), "Macrophage", "Tumour") }
#' @importFrom dbscan frNN
#' @export
plotProximity <- function(object, spe, genes, focal, neighbour,
                          radius = 30, breaks = c(0, 2, 4, 6, 8, Inf),
                          box_colour = "#4F8B5E", assay_name = "counts") {
  stopifnot(methods::is(object, "PACEFit"))
  ct     <- as.character(SummarizedExperiment::colData(spe)[[object@params$celltype_col]])
  coords <- as.matrix(SpatialExperiment::spatialCoords(spe))
  focal_idx <- which(ct == focal)
  nbr_idx   <- which(ct == neighbour)

  ## number of neighbour-type cells within `radius` um of each focal cell
  fr <- dbscan::frNN(x = coords[nbr_idx, , drop = FALSE], eps = radius,
                     query = coords[focal_idx, , drop = FALSE])
  n_neighbours <- lengths(fr$id)

  counts   <- SummarizedExperiment::assay(spe, assay_name)
  expr_mat <- as.matrix(counts[genes, focal_idx, drop = FALSE])   # genes x cells (RAW)
  plot_data <- data.frame(
    n_neighbours = rep(n_neighbours, times = length(genes)),
    gene = factor(rep(genes, each = length(focal_idx)), levels = genes),
    expr = as.numeric(t(expr_mat)))
  plot_data$bin <- cut(plot_data$n_neighbours, breaks = breaks, include.lowest = TRUE)
  plot_data <- plot_data[!is.na(plot_data$bin), ]
  ## relabel the open top bin with the observed maximum count
  lev <- levels(plot_data$bin)
  lev[length(lev)] <- sprintf("(%g,%d]", breaks[length(breaks) - 1L],
                              max(plot_data$n_neighbours))
  levels(plot_data$bin) <- lev

  p <- ggplot(plot_data, aes(.data$bin, .data$expr)) +
    geom_boxplot(colour = box_colour, fill = "white", linewidth = 0.6,
                 outlier.colour = box_colour, outlier.alpha = 0.25, outlier.size = 1.2) +
    stat_summary(fun = mean, geom = "point", colour = box_colour, size = 3.2) +
    labs(x = paste("Number of", neighbour, "neighbours"),
         y = paste(focal, "RAW counts")) +
    theme_minimal(base_size = 15) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank())

  if (length(genes) == 1L)
    p + labs(title = paste(genes, "expression"),
             y = paste(focal, genes, "counts"))
  else
    p + facet_wrap(~ gene, scales = "free_y")
}
