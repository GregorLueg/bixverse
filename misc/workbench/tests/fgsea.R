seed <- 10101

set.seed(seed)

gene_sets_no <- 50

ranks <- rnorm(100)
names(ranks) <- paste0("a", 1:100)

universe <- names(ranks)

gene_sets <- purrr::map(
  1:gene_sets_no,
  ~ {
    set.seed(seed * .x + 1)
    size <- sample(5:20, 1)
    sample(universe, size, replace = FALSE)
  }
)

names(gene_sets) <- purrr::map_chr(
  1:gene_sets_no,
  ~ {
    set.seed(seed + .x + 1)
    paste(sample(LETTERS, 3), collapse = "")
  }
)

any(duplicated(names(gene_sets)))

new_order <- names(sort(purrr::map_dbl(gene_sets, length)))
gene_sets <- gene_sets[new_order]

ranks_sorted <- sort(ranks)

ranks_sorted

gene_set_i <- gene_sets[[1]]

gene_set_i <- c("a93", "a46", "a33", "a27", "a35")

N_r <- sum(abs(ranks_sorted[gene_set_i]))

P_hit <- sum(ranks_sorted[gene_set_i]) / N_r

P_miss <- 1 / (100 - length(gene_set_i))

ES <- P_hit - P_miss
