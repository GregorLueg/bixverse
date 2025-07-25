# synthetic data ---------------------------------------------------------------

## gene expression version -----------------------------------------------------

### version 1 ------------------------------------------------------------------

#' Generates a synthetic, pseudo gene expression matrix
#'
#' @description
#' This is a function to create a synthetic, pseudo gene expression matrix with
#' a pre-defined number of total genes, groups and number of genes
#' differentially expressed in the groups. The function will create in a
#' standard setting a 1000 x 90 matrix, and add optionally a small group into
#' the mixture that is slightly more difficult to detect
#'
#' @param no_grps Integer. Number of initial larger groups. Default: 3L.
#' @param per_group Integer. Number of samples per larger groups. Default: 30L.
#' @param total_genes Integer. Number of total genes in the matrix. Default:
#' 1000L.
#' @param no_genes_up Integer. Number of genes per group that are differentially
#' expressed. Default: 100L.
#' @param add_small_group Boolean. Add a smaller group that overlaps with group
#' 1? Default: TRUE.
#' @param size_small_grp Integer. Size of the smaller group.
#' @param seed Integer. Initial random seed for generation of the synthetic
#' data. Default: 10101L.
#'
#' @return A `synthetic_matrix_simple` class containing:
#' \itemize{
#'  \item mat - The random matrix
#'  \item diff - List of differentially expressed genes per group
#'  \item group - Look-up vector for sample to group mapping
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' synthetic_GEX <- create_synthetic_signal_matrix()
#' synthetic_signal_mat <- synthetic_GEX$mat
#' }
synthetic_signal_matrix <- function(
  no_grps = 3,
  per_group = 30,
  total_genes = 1000L,
  no_genes_up = 100L,
  add_small_group = TRUE,
  size_small_grp = 5L,
  seed = 10101L
) {
  # hard coded
  mean_initial <- 5
  variance_initial <- 5
  mean_up <- 10
  variance_up <- 5

  # checks
  checkmate::qassert(no_grps, "R1(0,)")
  checkmate::qassert(per_group, "R1(0,)")
  checkmate::qassert(total_genes, "I1(0,)")
  checkmate::qassert(no_genes_up, "I1(0,)")
  checkmate::qassert(add_small_group, "B1")
  checkmate::qassert(size_small_grp, "I1(0,)")
  checkmate::qassert(seed, "I1")

  # matrix generation
  set.seed(seed)
  mat <- matrix(
    rnorm(
      total_genes * per_group * no_grps,
      mean_initial,
      variance_initial
    ),
    nrow = total_genes,
    ncol = per_group * no_grps
  )

  # add row and column names
  rownames(mat) <- paste0("gene", seq_len(total_genes))
  colnames(mat) <- paste0("sample", seq_len(per_group * no_grps))

  # create group factor
  group <- factor(unlist(lapply(seq_len(no_grps), function(x) {
    rep(paste0("group", x), per_group)
  })))
  names(group) <- colnames(mat)

  # create differentially expressed genes
  set.seed(seed + 1)
  diff <- lapply(seq_len(no_grps), function(x) {
    start_idx <- ((x - 1) * no_genes_up) + 1
    end_idx <- ((x - 1) * no_genes_up) + no_genes_up
    diff_genes <- rownames(mat)[start_idx:end_idx]

    # update matrix with differential expression
    group_mask <- group == paste0("group", x)
    mat[diff_genes, group_mask] <<-
      mat[diff_genes, group_mask] +
      rnorm(no_genes_up, mean_up, variance_up)

    return(diff_genes)
  })
  names(diff) <- paste0("group", seq_len(no_grps))

  # add small group if requested
  if (add_small_group) {
    small_group_genes <- (total_genes - no_genes_up):total_genes
    small_group_samples <- seq_len(size_small_grp)

    mat[small_group_genes, small_group_samples] <-
      mat[small_group_genes, small_group_samples] + mean_up

    diff$small_group <- rownames(mat)[small_group_genes]

    group <- as.character(group)
    group[small_group_samples] <- "small_group"
    names(group) <- colnames(mat)
    group <- factor(group)
  }

  # ensure no negative expression values
  mat[mat < 0] <- 0

  # finalise results
  result <- list(
    mat = mat,
    diff = diff,
    group = group
  )

  class(result) <- "synthetic_matrix_simple"

  return(result)
}

## contrastive pca synthetic data ----------------------------------------------

#' @title
#' Generates synthetic data for contrastive PCA exploration.
#'
#' @description
#' Generates three elements: target matrix, background matrix and labels for the
#' target matrix.
#'
#' @param seed Integer. Initial random seed for generation of the synthetic
#' data. Default: 10101L.
#'
#' @return A `cpca_synthetic_data` class with the following elements:
#' \itemize{
#'  \item target - The target matrix.
#'  \item background - The background matrix.
#'  \item target_labels - The target labels
#' }
#'
#' @importFrom magrittr `%>%`
#'
#' @export
synthetic_c_pca_data <- function(seed = 10101L) {
  # Checks
  checkmate::qassert(seed, "I1")
  set.seed(seed)

  # Background matrix
  background <- matrix(0, nrow = 400, ncol = 30)
  background[, 1:10] <- matrix(rnorm(400 * 10, mean = 0, sd = 10), ncol = 10)
  background[, 11:20] <- matrix(rnorm(400 * 10, mean = 0, sd = 3), ncol = 10)
  background[, 21:30] <- matrix(rnorm(400 * 10, mean = 0, sd = 1), ncol = 10)

  # Target matrix with 4 clusters
  target <- matrix(0, nrow = 400, ncol = 30)
  target[, 1:10] <- matrix(rnorm(400 * 10, mean = 0, sd = 10), ncol = 10)

  # Group 1
  target[1:100, 11:20] <- matrix(rnorm(100 * 10, mean = 0, sd = 1), ncol = 10)
  target[1:100, 21:30] <- matrix(rnorm(100 * 10, mean = 0, sd = 1), ncol = 10)

  # Group 2
  target[101:200, 11:20] <- matrix(rnorm(100 * 10, mean = 0, sd = 1), ncol = 10)
  target[101:200, 21:30] <- matrix(rnorm(100 * 10, mean = 3, sd = 1), ncol = 10)

  # Group 3
  target[201:300, 11:20] <- matrix(rnorm(100 * 10, mean = 6, sd = 1), ncol = 10)
  target[201:300, 21:30] <- matrix(rnorm(100 * 10, mean = 0, sd = 1), ncol = 10)

  # Group 4
  target[301:400, 11:20] <- matrix(rnorm(100 * 10, mean = 6, sd = 1), ncol = 10)
  target[301:400, 21:30] <- matrix(rnorm(100 * 10, mean = 3, sd = 1), ncol = 10)

  # Labels
  target_labels <- factor(as.character(rep(0:3, each = 100)))

  # Assign row and column names
  rownames(background) <- rownames(target) <- paste0("sample_", 1:400)
  colnames(background) <- colnames(target) <- paste0("feature_", 1:30)

  result <- list(
    target = t(target),
    background = t(background),
    target_labels = target_labels
  )

  class(result) <- "cpca_synthetic_data"

  return(result)
}

## gene module data ------------------------------------------------------------

#' @title
#' Generates synthetic gene module data.
#'
#' @description
#' Generates an artifical matrix of active modules (for a maximum of 8) in
#' samples being active. Also allows for overlap between the modules. Designed
#' to test some of the algorithms.
#'
#' @param n_samples Integer. Number of samples.
#' @param n_genes Integer. Number of genes.
#' @param n_modules Integer. Number of active gene modules. To a maximum of 10.
#' @param overlap_fraction Float. How much the same modules can be active
#' in the same samples.
#' @param seed Integer. Initial random seed for generation of the synthetic
#' data. Default: 10101L.
#'
#' @return A `synthetic_matrix_modules` class with the following items:
#' \itemize{
#'  \item data - The data matrix.
#'  \item metadata - The sample metadata.
#' }
#'
#' @importFrom magrittr `%>%`
#'
#' @export
generate_gene_module_data <- function(
  n_samples = 24L,
  n_genes = 60L,
  n_modules = 4L,
  overlap_fraction = 0.25,
  seed = 10101L
) {
  # checks
  checkmate::qassert(n_samples, "I1")
  checkmate::qassert(n_genes, "I1")
  checkmate::qassert(n_modules, "I1[1, 11)")
  checkmate::qassert(overlap_fraction, "N1(0, 1)")
  checkmate::qassert(seed, "I1")

  # body
  set.seed(seed)

  data <- matrix(0, nrow = n_samples, ncol = n_genes)

  genes_per_module <- n_genes %/% n_modules
  gene_boundaries <- seq(1, n_genes + 1, by = genes_per_module)
  if (length(gene_boundaries) > n_modules + 1) {
    gene_boundaries <- gene_boundaries[1:(n_modules + 1)]
  }
  gene_boundaries[n_modules + 1] <- n_genes + 1 # Ensure last boundary covers all genes

  cells_per_module <- n_samples %/% n_modules
  overlap_size <- round(cells_per_module * overlap_fraction)

  modules <- vector("list", n_modules)
  active_cells <- vector("list", n_modules)

  for (i in 1:n_modules) {
    gene_start <- gene_boundaries[i]
    gene_end <- gene_boundaries[i + 1] - 1

    modules[[i]] <- list(
      start = gene_start,
      end = gene_end
    )

    cell_start <- max(1, (i - 1) * cells_per_module - overlap_size + 1)
    cell_end <- min(n_samples, i * cells_per_module + overlap_size)

    active_cells[[i]] <- cell_start:cell_end
  }

  # Simulate module activities
  for (i in seq_along(modules)) {
    module <- modules[[i]]

    for (cell in active_cells[[i]]) {
      base_activity <- 1.0 + runif(1, -0.2, 0.2)

      for (gene in module$start:module$end) {
        data[cell, gene] <- base_activity + runif(1, -0.1, 0.1)
      }
    }
  }

  # Add background noise
  noise <- matrix(
    runif(n_samples * n_genes, -0.05, 0.05),
    nrow = n_samples,
    ncol = n_genes
  )
  data <- data + noise

  rownames(data) <- sprintf("sample_%i", 1:n_samples)
  colnames(data) <- sprintf("feature_%i", 1:n_genes)

  activity <- purrr::map(active_cells, \(sample_idx) {
    1:n_samples %in% sample_idx
  })

  activity_dt <- data.table::data.table(
    sample_id = sprintf("sample_%i", 1:n_samples)
  )

  for (i in seq_along(activity)) {
    activity_dt[, c(sprintf("module_%i_active", i)) := activity[[i]]]
  }

  result <- list(data = data, meta_data = activity_dt)

  class(result) <- "synthetic_matrix_modules"

  return(result)
}

# plots ------------------------------------------------------------------------

## simple data -----------------------------------------------------------------

#' Plot the (simple) synthetic gene expression
#'
#' @param synthetic_gex `synthetic_matrix_simple` class. Output from
#' [bixverse::synthetic_signal_matrix()].
#' @param ... Additional params
#'
#' @return A plotted heatmap showing the DEG.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' synthetic_gex <- synthetic_signal_matrix()
#'
#' plot(synthetic_gex)
#' }
plot.synthetic_matrix_simple <- function(synthetic_gex, ...) {
  # checks
  checkmate::assertClass(synthetic_gex, "synthetic_matrix_simple")

  # pull out relevant parts
  mat <- synthetic_gex$mat
  group <- synthetic_gex$group

  # plot
  heatmap(
    mat,
    Rowv = NA,
    Colv = NA,
    col = colorRampPalette(c("blue", "white", "red"))(100),
    scale = "none",
    ColSideColors = rainbow(length(levels(group)))[group],
    labCol = FALSE,
    labRow = FALSE,
    xlab = "Samples",
    ylab = "Features",
    main = "Synthetic gene expression data"
  )
}

## contrastive pca -------------------------------------------------------------

#' @title
#' Plot the contrastive PCA example data
#'
#' @param cpca_data `cpca_synthetic_data` class. Output from
#' [bixverse::synthetic_c_pca_data()].
#' @param ... Additional params
#'
#' @return A ggplot showing the two heatmaps from the target and background
#' matrix.
#'
#' @export
#'
#' @import patchwork
#' @import ggplot2
#' @importFrom magrittr `%>%`
plot.cpca_synthetic_data <- function(cpca_data, ...) {
  # checks
  checkmate::assertClass(cpca_data, "cpca_synthetic_data")

  # get the data
  target_df <- as.data.table(t(cpca_data$target), keep.rownames = "samples") %>%
    .[, grp := cpca_data$target_labels] %>%
    melt(id.vars = c("samples", "grp"))

  background_df <- as.data.table(
    t(cpca_data$background),
    keep.rownames = "samples"
  ) %>%
    melt(id.vars = "samples")

  # plots
  p1 <- ggplot(
    data = target_df,
    mapping = aes(y = samples, x = variable, fill = value)
  ) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    ggtitle("Target data") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    ) +
    xlab("Features") +
    ylab("Samples") +
    geom_rug(
      data = unique(target_df[, c("samples", "grp")]),
      mapping = aes(y = samples, color = grp),
      inherit.aes = FALSE,
      sides = "l",
      linewidth = 2,
      length = unit(0.0125, "npc")
    ) +
    labs(fill = "Value:", colour = "Group:")

  p2 <- ggplot(
    data = background_df,
    mapping = aes(y = samples, x = variable, fill = value)
  ) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    ggtitle("Background data") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    ) +
    xlab("Features") +
    ylab("Samples") +
    labs(fill = "Value:")

  p_final <- p1 +
    p2 +
    patchwork::plot_annotation(
      title = "Contrastive PCA",
      subtitle = "Synthetic data"
    )

  p_final
}
