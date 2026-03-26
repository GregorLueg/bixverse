# mitch test -------------------------------------------------------------------

## data ------------------------------------------------------------------------

set.seed(42L)

contrast_data <- matrix(data = rnorm(3 * 26), nrow = 26)
colnames(contrast_data) <- sprintf("contrast_%i", 1:3)
rownames(contrast_data) <- letters

gene_sets <- purrr::map(
  1:4,
  ~ {
    sample(x = letters, size = .x + 3)
  }
)
names(gene_sets) <- sprintf("pathway_%s", LETTERS[1:4])

## tests -----------------------------------------------------------------------

res <- calc_mitch(
  contrast_mat = contrast_data,
  gene_set_list = gene_sets
)

expect_equal(
  current = res$pathway_names,
  target = c("pathway_C", "pathway_D", "pathway_B"),
  info = "mitch - expected pathway names"
)

expect_equal(
  current = res$pathway_sizes,
  target = c(6, 7, 5),
  info = "mitch - expected pathway sizes"
)

expect_equal(
  current = res$manova_pval,
  target = c(0.006953769, 0.486574059, 0.660185031),
  eps = 1e-7,
  info = "mitch - expected manova pvals"
)

### with NA situation ----------------------------------------------------------

x_coords <- sample(nrow(contrast_data), size = 4, replace = TRUE)
y_coords <- sample(ncol(contrast_data), size = 4, replace = TRUE)

bad_contrast_data <- contrast_data

for (i in seq_along(x_coords)) {
  bad_contrast_data[x_coords[i], y_coords[i]] <- NA
}

expect_error(
  current = calc_mitch(
    contrast_mat = bad_contrast_data,
    gene_set_list = gene_sets
  ),
  info = paste("mitch - throw error with NAs")
)

### special case - 2 contrasts -------------------------------------------------

contrast_data_2 <- contrast_data[, 1:2]

res_2 <- calc_mitch(
  contrast_mat = contrast_data_2,
  gene_set_list = gene_sets
)

expect_equal(
  current = res_2$pathway_names,
  target = c("pathway_C", "pathway_D", "pathway_B"),
  info = "mitch - expected pathway names (only two contrasts)"
)

expect_equal(
  current = res_2$pathway_sizes,
  target = c(6, 7, 5),
  info = "mitch - expected pathway sizes (only two contrasts)"
)

expect_equal(
  current = res_2$manova_pval,
  target = c(0.03262419, 0.37960291, 0.75588214),
  eps = 1e-7,
  info = "mitch - expected manova pvals (only two contrasts)"
)

## direct comparison mitch -----------------------------------------------------

local({
  if (requireNamespace("mitch", quietly = TRUE)) {
    # load in the mitch example data
    data(myImportedData, genesetsExample, package = "mitch")

    mitch_res <- suppressMessages(mitch::mitch_calc(
      myImportedData,
      genesetsExample,
      priority = 'significance',
      minsetsize = 5,
      cores = 2
    ))

    bixverse_res <- calc_mitch(
      contrast_mat = as.matrix(myImportedData),
      gene_set_list = genesetsExample
    )

    expect_equal(
      current = bixverse_res$pathway_names,
      target = mitch_res$enrichment_result$set,
      info = "mitch direct comparison - pathways"
    )

    expect_equal(
      current = bixverse_res$manova_fdr,
      target = mitch_res$enrichment_result$p.adjustMANOVA,
      info = "mitch direct comparison - adjusted MANOVA pvalues"
    )

    expect_equal(
      current = bixverse_res$s.rna,
      target = mitch_res$enrichment_result$s.rna,
      info = "mitch direct comparison - s.rna"
    )

    expect_equal(
      current = bixverse_res$s.k9a,
      target = mitch_res$enrichment_result$s.k9a,
      info = "mitch direct comparison - s.k9a"
    )

    expect_equal(
      current = bixverse_res$s.k36a,
      target = mitch_res$enrichment_result$s.k36a,
      info = "mitch direct comparison - s.k36a"
    )

    expect_equal(
      current = bixverse_res$p.rna,
      target = mitch_res$enrichment_result$p.rna,
      info = "mitch direct comparison - p.rna"
    )

    expect_equal(
      current = bixverse_res$p.k9a,
      target = mitch_res$enrichment_result$p.k9a,
      info = "mitch direct comparison - p.k9a"
    )

    expect_equal(
      current = bixverse_res$p.k36a,
      target = mitch_res$enrichment_result$p.k36a,
      info = "mitch direct comparison - p.k36a"
    )
  }
})
