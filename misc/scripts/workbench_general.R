n_samples = 24L
n_genes = 60L
n_modules = 4L
overlap_fraction = 0.25
seed = 10101L

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

results <- list(data = data, meta_data = activity_dt)
