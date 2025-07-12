go_data <- load_go_human_data()

go_genes <- go_data$go_to_genes

go_genes_ls <- split(go_genes$ensembl_id, go_genes$go_id)

go_genes_ls <- purrr::keep(go_genes_ls, \(x) length(x) > 3L)

length(go_genes_ls)

results <- gse_hypergeometric(
  target_genes = target_genes_1,
  gene_set_list = go_genes_ls
)


go_parent_child_dt <- go_data$gene_ontology[
  relationship %in% c("is_a", "part_of")
] %>%
  setnames(
    old = c("from", "to", "relationship"),
    new = c("parent", "child", "type")
  )

results_simplified <- simplify_hypergeom_res(
  res = results,
  parent_child_dt = go_parent_child_dt,
  weights = setNames(c(0.8, 0.6), c("is_a", "part_of"))
)
