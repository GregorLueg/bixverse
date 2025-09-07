# library(mitch)

# ?mitch_calc

# data(myImportedData, genesetsExample)

# tictoc::tic()
# resExample <- mitch_calc(
#   myImportedData,
#   genesetsExample,
#   priority = 'significance',
#   minsetsize = 5,
#   cores = 2
# )
# tictoc::toc()

# microbenchmark::microbenchmark(
#   r = {
#     mitch_calc(
#       myImportedData,
#       genesetsExample,
#       priority = 'significance',
#       minsetsize = 5,
#       cores = 2
#     )
#   },
#   rust = {
#     rs_mitch_calc(
#       as.matrix(myImportedData),
#       rownames(myImportedData),
#       genesetsExample,
#       5
#     )
#   },
#   times = 50L
# )

# tictoc::tic()
# rs_result <- rs_mitch_calc(
#   as.matrix(myImportedData),
#   rownames(myImportedData),
#   genesetsExample,
#   5
# )
# tictoc::toc()

# rs_result$pathway_names[order(rs_result$manova_pvals)]
# rs_result$pathway_sizes[order(rs_result$manova_pvals)]
# rs_result$sd[order(rs_result$manova_pvals)]

# resExample$enrichment_result

# input_profile <- myImportedData
# genesets = genesetsExample
# minsetsize = 5
# cores = 2
# priority = 'effect'

# ranked_profile <- mitch:::mitch_rank(input_profile)

# enrichment_result <- mitch:::MANOVA(
#   ranked_profile,
#   genesets,
#   minsetsize <- minsetsize,
#   cores = cores,
#   priority = priority
# )

# x = ranked_profile

# sets <- names(genesets)

# hypotenuse <- function(x) {
#   sqrt(sum(
#     unlist(lapply(x, function(x) {
#       x^2
#     })),
#     na.rm = TRUE
#   ))
# }

# x <- scord

# sqrt(sum(x^2))

# lapply(x, function(x) {
#   x^2
# })

# HYPOT <- hypotenuse(apply(x, 2, length))

# set = "Biological oxidations"

# inset <- rownames(x) %in% as.character(unlist(genesets[set]))

# rextendr::document()

# rs_fast_mitch(x, inset)

# fit <- manova(x ~ inset)
# sumMANOVA <- summary.manova(fit)
# sumAOV <- summary.aov(fit)
# pMANOVA <- sumMANOVA$stats[1, "Pr(>F)"]

# summary(fit)$stats

# NROW <- nrow(x)

# raov <- lapply(sumAOV, function(zz) {
#   zz[1, "Pr(>F)"]
# })
# raov <- unlist(raov)
# names(raov) <- gsub("^ Response ", "p.", names(raov))
# # S coordinates
# NOTINSET <- colMeans(x[!inset, ], na.rm = TRUE)
# scord <- (2 * (colMeans(x[inset, ], na.rm = TRUE) - NOTINSET)) / NROW
# names(scord) <- paste0("s-", names(scord))
# # calculate the hypotenuse length of s scores
# s.dist <- hypotenuse(scord)
# names(s.dist) <- "s.dist"
# mysd <- sd(scord)
# names(mysd) <- "SD"

# data.frame(
#   set,
#   setSize = sum(inset),
#   pMANOVA,
#   t(scord),
#   t(raov),
#   t(s.dist),
#   t(mysd),
#   stringsAsFactors = FALSE
# )
