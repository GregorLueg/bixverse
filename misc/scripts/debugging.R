original_function <- function(
  stats,
  selectedStats,
  gseaParam = 1,
  returnAllExtremes = FALSE,
  returnLeadingEdge = FALSE,
  scoreType = c("std", "pos", "neg")
) {
  scoreType <- match.arg(scoreType)
  S <- selectedStats
  r <- stats
  p <- gseaParam
  S <- sort(S)
  stopifnot(all(head(S, -1) < tail(S, -1)))
  m <- length(S)
  N <- length(r)
  if (m == N) {
    stop("GSEA statistic is not defined when all genes are selected")
  }
  NR <- (sum(abs(r[S])^p))
  rAdj <- abs(r[S])^p
  if (NR == 0) {
    rCumSum <- seq_along(rAdj) / length(rAdj)
  } else {
    rCumSum <- cumsum(rAdj) / NR
  }
  tops <- rCumSum - (S - seq_along(S)) / (N - m)
  print("tops:")
  print(tops)
  if (NR == 0) {
    bottoms <- tops - 1 / m
  } else {
    bottoms <- tops - rAdj / NR
  }
  maxP <- max(tops)
  minP <- min(bottoms)
  switch(
    scoreType,
    std = geneSetStatistic <- ifelse(
      maxP == -minP,
      0,
      ifelse(maxP > -minP, maxP, minP)
    ),
    pos = geneSetStatistic <- maxP,
    neg = geneSetStatistic <- minP
  )
  if (!returnAllExtremes && !returnLeadingEdge) {
    return(geneSetStatistic)
  }
  res <- list(res = geneSetStatistic)
  if (returnAllExtremes) {
    res <- c(res, list(tops = tops, bottoms = bottoms))
  }
  if (returnLeadingEdge) {
    switch(
      scoreType,
      std = leadingEdge <- if (maxP > -minP) {
        S[seq_along(S) <= which.max(tops)]
      } else if (maxP < -minP) {
        rev(S[seq_along(S) >= which.min(bottoms)])
      } else {
        NULL
      },
      pos = leadingEdge <- S[seq_along(S) <= which.max(tops)],
      neg = leadingEdge <- rev(S[seq_along(S) >= which.min(bottoms)])
    )
    res <- c(res, list(leadingEdge = leadingEdge))
  }
  res
}

fgsea_res <- original_function(
  stats = stats,
  selectedStats = pathway_indices_r$pathway_pos,
  returnLeadingEdge = TRUE,
  returnAllExtremes = TRUE
)


fgsea_res$bottoms
