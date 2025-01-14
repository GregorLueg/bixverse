# Load in and extend the Rust documentation...
rextendr::document()
devtools::document()
devtools::load_all()

devtools::install()

library(magrittr)

protein_coding_genes = data.table::fread("~/Desktop/protein_coding_genes.csv")

seed = 123
set.seed(seed)

universe = protein_coding_genes$id

gene_sets = purrr::map(1:5000, ~{
  set.seed(seed + .x + 1)
  size = sample(20:100, 1)
  sample(universe, size, replace = FALSE)
})

target_gene_sets = purrr::map(1:1000, ~{
  set.seed(.x * seed)
  size = sample(50:100, 1)
  sample(universe, size, replace = FALSE)
})

target_genes = target_gene_sets[[2]]

tictoc::tic()
test = hypergeom_test(
  target_genes = target_genes,
  gene_sets = gene_sets,
  gene_universe = universe
)
tictoc::toc()


tictoc::tic()
test_2 = hypergeom_test_list(
  target_gene_lists = target_gene_sets,
  gene_sets = gene_sets,
  gene_universe = universe
)
tictoc::toc()


future::plan(future::multisession(workers = 5))

tictoc::tic()
test_3 = furrr::future_map(
  target_gene_sets,
  hypergeom_test,
  gene_sets = gene_sets,
  gene_universe = universe,
  .progress = T
)
tictoc::toc()

future::plan(future::sequential())
gc()

do.call(cbind, test_3[[1]])



# Compare to Aether ----

# Old Aether code...
.gse_hypergeom_test <- function(target_genes,
                                gene_set,
                                gene_universe,
                                return_contingency_info = F) {
  # Checks
  checkmate::qassert(target_genes, "S+")
  checkmate::qassert(gene_set, "S+")
  checkmate::qassert(gene_universe, "S+")
  checkmate::qassert(return_contingency_info, "B1")
  # Function body
  gene_universe_length <- length(gene_universe)
  gene_set_length <- length(intersect(gene_set, gene_universe))
  trials <- length(target_genes)
  hits <- length(intersect(target_genes, gene_set))

  # Bizarrely enough, it's faster to do this like this instead of using fisher.test?
  pval <- phyper(
    hits - 1,
    gene_set_length,
    gene_universe_length - gene_set_length,
    trials,
    lower.tail = F
  )
  # Contingency table A = in target genes; B = in gene set
  # In target genes and in gene set
  A1_B1 <- hits
  # Not in target genes but in gene set
  A0_B1 <- gene_set_length - hits
  # In target genes but not in gene set
  A1_B0 <- trials - hits
  # Not present in either
  A0_B0 <- gene_universe_length - gene_set_length - trials + hits

  OR <- (A1_B1 / A1_B0) / (A0_B1 / A0_B0)

  if (return_contingency_info) {
    return(c(
      pval = pval,
      OR = OR,
      A1_B1 = A1_B1,
      A0_B1 = A0_B1,
      A1_B0 = A1_B0,
      A0_B0 = A0_B0
    ))
  } else {
    return(c(pval = pval, OR = OR))
  }
}

GSE_hypergeometric <- function(target_genes,
                               gene_set_list,
                               gene_universe = NULL,
                               threshold = 0.05,
                               return_contingency_info = F,
                               verbose = T) {
  # Checks
  checkmate::qassert(target_genes, "S+")
  checkmate::assertList(gene_set_list, types = "character")
  checkmate::qassert(gene_universe, c("0", "S+"))
  checkmate::qassert(verbose, "B1")
  checkmate::qassert(return_contingency_info, "B1")
  checkmate::qassert(threshold, c("R1[0,1]", "0"))
  # Function body
  if (is.null(gene_universe)) {
    if (verbose) message("No gene universe given. Function will use the represented genes in the pathways/gene sets as reference.")
    gene_universe <- unique(unlist(gene_set_list))
  }
  # Throw warnings if there are duplicates in the target genes
  if(any(duplicated(target_genes))) warning("Duplicated genes found in target genes.")
  gse_results <- data.table::data.table(
    name = names(gene_set_list),
    purrr::map_dfr(
      gene_set_list,
      .gse_hypergeom_test,
      target_genes = target_genes,
      gene_universe = gene_universe,
      return_contingency_info = return_contingency_info
    )
  ) %>%
    .[, FDR := p.adjust(pval, method = "BH")]

  # Additional filter
  if(!is.null(threshold)) gse_results = gse_results[FDR <= threshold]

  # Add additional columns
  if(return_contingency_info) {
    gse_results[, `:=`(query = A1_B1 + A1_B0,
                       reference = A1_B1 + A0_B1,
                       intersection = A1_B1)]
  }

  gse_results
}

GSE_hypergeometric_over_list <- function(target_gene_list,
                                         gene_set_list,
                                         gene_universe = NULL,
                                         threshold = 0.05,
                                         verbose = T,
                                         parallel = F,
                                         return_contingency_info = F,
                                         return_as_DT = T,
                                         no_cores = as.integer(parallel::detectCores() / 2)) {
  # Checks
  checkmate::assertList(target_gene_list, types = "character")
  checkmate::assertList(gene_set_list, types = "character")
  checkmate::qassert(gene_universe, c("0", "S+"))
  checkmate::qassert(threshold, "R1[0,1]")
  checkmate::qassert(no_cores, "I1")
  checkmate::qassert(parallel, "B1")
  checkmate::qassert(return_contingency_info, "B1")
  checkmate::qassert(return_as_DT, "B1")
  checkmate::qassert(verbose, "B1")
  # Function body
  if (is.null(gene_universe)) {
    if (verbose)
      message(
        "No gene universe given. Function will use the represented genes in the pathways/gene sets as reference."
      )
    gene_universe <- unique(unlist(gene_set_list))
  }
  if (!parallel) {
    if (verbose) {
      message("Using sequential calculations.")
    }
    results <- foreach(i = seq_along(target_gene_list)) %do% {
      enrichment_result_i <- GSE_hypergeometric(target_gene_list[[i]],
                                                gene_set_list,
                                                gene_universe = gene_universe,
                                                threshold = threshold,
                                                return_contingency_info = return_contingency_info)
    }
  } else {
    if (verbose) {
      message("Using parallel calculations with ", no_cores, " cores.")
    }
    # I will keep the doSNOW back-end for the time being, even knowing it causes issues in times
    # Progress bar
    cl <- parallel::makeCluster(no_cores)
    registerDoSNOW(cl)
    pb <- txtProgressBar(max = length(target_gene_list),
                         style = 3)
    progress <- function(n) {
      setTxtProgressBar(pb, n)
    }
    opts <- list(progress = progress)
    # Iterate through everything
    results <- foreach(
      i = seq_along(target_gene_list),
      .packages = c("data.table"),
      .options.snow = opts
    ) %dopar% {
      enrichment_result_i <-
        Aether::GSE_hypergeometric(
          target_gene_list[[i]],
          gene_set_list,
          gene_universe = gene_universe,
          threshold = threshold,
          return_contingency_info = return_contingency_info
        )
    }
    parallel::stopCluster(cl)
  }
  on.exit(gc()) # Should clean up all connections if things break...
  names(results) <- names(target_gene_list)

  if(return_as_DT) results = purrr::imap_dfr(results, ~{ .x[, module_name := .y]})

  results
}

tictoc::tic()
GSE_hypergeometric(target_genes = target_genes,
                   gene_set_list = gene_sets)
tictoc::toc()

library(doSNOW)

tictoc::tic()
GSE_hypergeometric_over_list(
  target_gene_list = target_gene_sets,
  gene_set_list = gene_sets,
  parallel = TRUE
)
tictoc::toc()
