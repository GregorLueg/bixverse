# synthetic data ----

## general signal data ----

#' Generates a synthetic, pseudo gene expression matrix
#'
#' @description
#' This is a function to create a synthetic, pseudo gene expression matrix with
#' a pre-defined number of total genes, groups and number of genes
#' differentially expressed in the groups. The function will create in a
#' standard setting a 1000 x 90 matrix, and add a small group into the mixture
#' that is more difficult to detect.
#'
#' @param no_grps Integer. Number of initial larger groups. Default: 3L.
#' @param per_group Integer. Number of samples per larger groups. Default: 30L.
#' @param total_genes Integer. Number of total genes in the matrix.
#' Default: 1000L.
#' @param mean_initial Float. Initial mean expression. Default: 5.
#' @param variance_initial Float. Initial variance in the expression.
#' Default: 5.
#' @param mean_up Float. Mean expression of the genes that are differentially
#' expressed. Default: 10.
#' @param variance_up Float. Variance of the genes that are differentially
#' expressed. Default: 5.
#' @param no_genes_up Integer. Number of genes per group that are differentially
#' expressed. Default: 100L.
#' @param add_small_group Boolean. Add a smaller group that overlaps with group
#' 1? Default: TRUE.
#' @param size_small_grp Integer. Size of the smaller group.
#' @param seed Integer. Initial random seed for generation of the synthetic
#' data. Default: 10101L.
#'
#' @return A list containing:
#'   \itemize{
#'     \item mat - The random matrix
#'     \item diff - List of differentially expressed genes per group
#'     \item group - Look-up vector for sample to group mapping
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' synthetic_GEX <- create_synthetic_signal_matrix()
#' synthetic_signal_mat <- synthetic_GEX$mat
#' }
create_synthetic_signal_matrix <- function(no_grps = 3,
                                           per_group = 30,
                                           total_genes = 1000L,
                                           mean_initial = 5,
                                           variance_initial = 5,
                                           mean_up = 10,
                                           variance_up = 5,
                                           no_genes_up = 100L,
                                           add_small_group = TRUE,
                                           size_small_grp = 5L,
                                           seed = 10101L) {
  # Input validation
  checkmate::qassert(no_grps, "R1(0,)")
  checkmate::qassert(per_group, "R1(0,)")
  checkmate::qassert(total_genes, "I1(0,)")
  checkmate::qassert(mean_initial, "R1")
  checkmate::qassert(variance_initial, "R1")
  checkmate::qassert(mean_up, "R1")
  checkmate::qassert(variance_up, "R1")
  checkmate::qassert(no_genes_up, "I1(0,)")
  checkmate::qassert(add_small_group, "B1")
  checkmate::qassert(size_small_grp, "I1(0,)")
  checkmate::qassert(seed, "I1")

  # Create initial random matrix
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

  # Add row and column names
  rownames(mat) <- paste0("gene", seq_len(total_genes))
  colnames(mat) <- paste0("sample", seq_len(per_group * no_grps))

  # Create group factor
  group <- factor(unlist(lapply(seq_len(no_grps), function(x) {
    rep(paste0("group", x), per_group)
  })))
  names(group) <- colnames(mat)

  # Create differentially expressed genes
  set.seed(seed + 1)
  diff <- lapply(seq_len(no_grps), function(x) {
    start_idx <- ((x - 1) * no_genes_up) + 1
    end_idx <- ((x - 1) * no_genes_up) + no_genes_up
    diff_genes <- rownames(mat)[start_idx:end_idx]

    # Update matrix with differential expression
    group_mask <- group == paste0("group", x)
    mat[diff_genes, group_mask] <<-
      mat[diff_genes, group_mask] +
      rnorm(no_genes_up, mean_up, variance_up)

    return(diff_genes)
  })
  names(diff) <- paste0("group", seq_len(no_grps))

  # Add small group if requested
  if (add_small_group) {
    small_group_genes <- (total_genes - no_genes_up):total_genes
    small_group_samples <- seq_len(size_small_grp)

    mat[small_group_genes, small_group_samples] <-
      mat[small_group_genes, small_group_samples] + mean_up

    diff$small_group <- rownames(mat)[small_group_genes]

    # Update group factor to include small group
    group <- as.character(group)
    group[small_group_samples] <- "small_group"
    names(group) <- colnames(mat)
    group <- factor(group)
  }

  # Ensure no negative expression values
  mat[mat < 0] <- 0

  # Prepare result
  result <- list(
    mat = mat,
    diff = diff,
    group = group
  )
  attr(result, "synthetic_example") <- TRUE

  return(result)
}

## cpca synthetic data ----

#' Generates synthetic data for cPCA exploration.
#'
#' @description Generates three elements: target matrix, background matrix
#' and labels for the target matrix.
#'
#' @param seed Integer. Initial random seed for generation of the synthetic data. Default: 10101L.
#'
#' @return A list. First elements is the target matrix, second element the background matrix and the
#' third one being a factor with the cluster labels of the target matrix.
#' @importFrom magrittr `%>%`
#'
#' @export
create_synthetic_cPCA_data <- function(seed = 10101L) {
  # Checks
  checkmate::qassert(seed, "I1")
  set.seed(seed)

  # Background matrix
  background <- matrix(0, nrow = 400, ncol = 30)
  background[, 1:10]  <- matrix(rnorm(400 * 10, mean = 0, sd = 10), ncol = 10)
  background[, 11:20] <- matrix(rnorm(400 * 10, mean = 0, sd = 3), ncol = 10)
  background[, 21:30] <- matrix(rnorm(400 * 10, mean = 0, sd = 1), ncol = 10)

  # Target matrix with 4 clusters
  target <- matrix(0, nrow = 400, ncol = 30)
  target[, 1:10]  <- matrix(rnorm(400 * 10, mean = 0, sd = 10), ncol = 10)

  # Group 1
  target[1:100, 11:20]  <- matrix(rnorm(100 * 10, mean = 0, sd = 1), ncol = 10)
  target[1:100, 21:30]  <- matrix(rnorm(100 * 10, mean = 0, sd = 1), ncol = 10)

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

  list(
    target = t(target),
    background = t(background),
    target_labels = target_labels
  )
}

# plotting functions ----

#' Plot the synthetic gene expression
#'
#' @param synthetic_GEX List. Output from Hephaestus::create_synthetic_signal_matrix().
#' It will throw an error if this is not an output from the function.
#'
#' @return A plotted heatmap showing the DEG.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' synthetic_GEX <- create_synthetic_signal_matrix()
#'
#' plot_synthetic_GEX_HT(synthetic_GEX)
#' }
plot_synthetic_GEX_HT <- function(synthetic_GEX) {
  # Checks
  # checkmate::assert(.checkSyntheticGEX(synthetic_GEX))
  # Get relevant parts
  mat <- synthetic_GEX$mat
  group <- synthetic_GEX$group
  # Plot
  heatmap(
    mat,
    Rowv = NA,
    Colv = NA,
    col = colorRampPalette(c("blue", "white", "red"))(100),
    scale = "none",
    ColSideColors = rainbow(length(levels(group)))[group],
    labCol = FALSE,
    labRow = FALSE
  )
}
