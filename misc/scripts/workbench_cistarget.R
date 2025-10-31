# test out the Rust version of cis targets -------------------------------------

rankins <- read_motif_ranking(
  "~/Downloads/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
)

geneList1 <- read.table(
  "~/repos/other/RcisTarget/inst/examples/hypoxiaGeneSet.txt"
)[, 1]

geneLists <- list(geneListName = geneList1)

geneLists_up <- list(
  geneListName = which(rownames(rankins) %in% geneLists$geneListName)
)

rextendr::document()

tictoc::tic()
res <- rs_cistarget(
  rankings = rankins,
  gs_list = geneLists_up,
  auc_threshold = as.integer(0.05 * nrow(rankins)),
  nes_threshold = 3.0,
  max_rank = nrow(rankins),
  method = "sth",
  n_mean = 10L
)
tictoc::toc()

res_dt <- data.table(
  motif = colnames(rankins)[res[[1]]$motif_idx],
  nes = res[[1]]$nes
)

setorder(res_dt, -nes)

res[[1]]$nes
