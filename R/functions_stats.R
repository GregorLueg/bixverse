# stat helpers -----------------------------------------------------------------

## topological overlap ---------------------------------------------------------

#' Calculate the TOM from a correlation matrix
#'
#' @param cor_mat Numerical matrix. The symmetric correlation matrix.
#' @param signed Boolean. Do you want to calculate the signed version. If set
#' to `FALSE`, the absolute correlation coefficients will be used.
#' @param version String. One of `c("v1", "v2")`. Defaults to `"v1"`.
#'
#' @details Calculates the topological overlap matrix from a correlation matrix.
#' The TOM is defined as:
#'
#' **Unsigned, v1:**
#'
#' \deqn{TOM_{ij} = \frac{a_{ij} + \sum_{k \neq i,j} a_{ik} a_{kj}}{\min(k_i, k_j) + 1 - a_{ij}}}
#'
#' **Signed, v1:**
#'
#' \deqn{TOM_{ij} = \frac{a_{ij} + \sum_{k \neq i,j} a_{ik} a_{kj}}{\min(k_i, k_j) + 1 - \left|a_{ij}\right|}}
#'
#' **Unsigned, v2:**
#'
#' \deqn{TOM_{ij} = 0.5 \left( a_{ij} + \frac{\sum_{k \neq i,j} a_{ik} a_{kj}}{\min(k_i, k_j) + a_{ij}} \right)}
#'
#' **Signed, v2:**
#'
#' \deqn{TOM_{ij} = 0.5 \left( a_{ij} + \frac{\sum_{k \neq i,j} a_{ik} a_{kj}}{\min(k_i, k_j) + \left|a_{ij}\right|} \right)}
#'
#' where \eqn{a_{ij}} is the affinity between nodes \eqn{i} and \eqn{j},
#' and \eqn{k_i = \sum_j a_{ij}} is the connectivity of node \eqn{i}.
#' For signed networks, connectivity is calculated as \eqn{k_i = \sum_j \left|a_{ij}\right|}.
#'
#' Version 2 uses a different normalisation approach that scales the shared
#' neighbour contribution separately before combining it with the direct
#' connection strength.
#'
#' @returns A symmetric matrix of the same dimensions as `cor_mat` containing
#' the topological overlap measures.
#'
#' @export
calculate_tom <- function(cor_mat, signed, version = c("v1", "v2")) {
  version <- match.arg(version)

  # checks
  checkmate::assertMatrix(cor_mat, nrows = ncol(cor_mat), ncols = nrow(cor_mat))
  checkmate::qassert(signed, "B1")
  checkmate::assertChoice(version, c("v1", "v2"))

  # body
  if (!signed) {
    cor_mat <- abs(cor_mat)
  }

  tom_mat <- rs_tom(x = cor_mat, tom_type = version, signed = signed)

  return(tom_mat)
}


#' Calculate the TOM from an expression matrix
#'
#' @param x Numerical matrix. The expression matrix. Assumes that columns are
#' the genes, and rows the samples.
#' @param signed Boolean. Do you want to calculate the signed version. If set
#' to `FALSE`, the absolute correlation coefficients will be used.
#' @param version String. One of `c("v1", "v2")`. Defaults to `"v1"`
#' @param cor_method String. One of `c("pearson", spearman)`. Defaults to
#' `"pearson"`.
#'
#' @details Calculates the topological overlap matrix from an expression matrix.
#' It will first calculate the specified correlation matrix and then generate
#' the TOM. The TOM is defined as:
#'
#' **Unsigned, v1:**
#'
#' \deqn{TOM_{ij} = \frac{a_{ij} + \sum_{k \neq i,j} a_{ik} a_{kj}}{\min(k_i, k_j) + 1 - a_{ij}}}
#'
#' **Signed, v1:**
#'
#' \deqn{TOM_{ij} = \frac{a_{ij} + \sum_{k \neq i,j} a_{ik} a_{kj}}{\min(k_i, k_j) + 1 - \left|a_{ij}\right|}}
#'
#' **Unsigned, v2:**
#'
#' \deqn{TOM_{ij} = 0.5 \left( a_{ij} + \frac{\sum_{k \neq i,j} a_{ik} a_{kj}}{\min(k_i, k_j) + a_{ij}} \right)}
#'
#' **Signed, v2:**
#'
#' \deqn{TOM_{ij} = 0.5 \left( a_{ij} + \frac{\sum_{k \neq i,j} a_{ik} a_{kj}}{\min(k_i, k_j) + \left|a_{ij}\right|} \right)}
#'
#' where \eqn{a_{ij}} is the affinity between nodes \eqn{i} and \eqn{j},
#' and \eqn{k_i = \sum_j a_{ij}} is the connectivity of node \eqn{i}.
#' For signed networks, connectivity is calculated as \eqn{k_i = \sum_j \left|a_{ij}\right|}.
#'
#' Version 2 uses a different normalisation approach that scales the shared
#' neighbour contribution separately before combining it with the direct
#' connection strength.
#'
#' @returns The topological overlap matrix.
#'
#' @export
calculate_tom_from_exp <- function(x, signed, version, cor_method) {
  # checks
  checkmate::assertMatrix(x)
  checkmate::qassert(signed, "B1")
  checkmate::assertChoice(version, c("v1", "v2"))

  spearman <- cor_method == "spearman"

  cor_mat <- rs_cor(x = x, spearman = spearman)

  # body
  if (!signed) {
    cor_mat <- abs(cor_mat)
  }

  tom_mat <- rs_tom(x = cor_mat, tom_type = version, signed = signed)

  return(tom_mat)
}

## F1 scores on confusion matrix -----------------------------------------------

#' F1 scores on top of a confusion matrix
#'
#' @description
#' Helper function to check for expected clustering vs actual clustering.
#'
#' @param clusters_a String or factor. The clustering of algorithm 1.
#' @param clusters_b String or factor. The clustering of algorithm 2.
#'
#' @returns Named vector with the F1 scores between the two clustering
#' algorithms.
#'
#' @export
f1_score_confusion_mat <- function(clusters_a, clusters_b) {
  # checks
  len_a <- length(clusters_a)
  len_b <- length(clusters_b)
  checkmate::qassert(clusters_a, sprintf("A%i", len_b))
  checkmate::qassert(clusters_b, sprintf("A%i", len_a))

  # function
  cm <- table(clusters_a, clusters_b)

  best_match <- apply(cm, 1, which.max)

  f1_scores <- sapply(seq_len(nrow(cm)), function(i) {
    tp <- cm[i, best_match[i]]
    fp <- sum(cm[, best_match[i]]) - tp
    fn <- sum(cm[i, ]) - tp

    precision <- tp / (tp + fp)
    recall <- tp / (tp + fn)

    2 * precision * recall / (precision + recall)
  })

  names(f1_scores) <- rownames(cm)

  return(f1_scores)
}
