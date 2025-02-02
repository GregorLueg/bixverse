library(magrittr)
library(data.table)

devtools::document()
rextendr::document()
devtools::load_all()
devtools::check()

# Community detections ----

edge_data = arrow::read_parquet("~/Desktop/biomind_downloads/processed_data/edges_OT_interactions.parquet") %>%
  as.data.table() %>%
  .[network_resource == 'STRING']

edge_data_clean = edge_data %>%
  .[interactionScore >= 0.85] %>%
  dplyr::select(from = `:START_ID`, to = `:END_ID`)

test_class = network_diffusions(edge_data_clean, weighted = FALSE, directed = TRUE)

print(test_class)

test_class

test_class@params

get_params(test_class)

set.seed(123)
genes = sample(igraph::V(test_class@graph)$name, 10)
genes.2 = sample(igraph::V(test_class@graph)$name, 10)
genes.3 = sample(igraph::V(test_class@graph)$name, 25)
diffusion_vector = rep(1, 10) %>% `names<-`(genes)
diffusion_vector.2 = rep(1, 10) %>% `names<-`(genes.2)

test_class <- diffuse_seed_nodes(test_class, diffusion_vector, 'max')

get_params(test_class, TRUE, TRUE)

test_class <- tied_diffusion(
  network_diffusions = test_class,
  diffusion_vector.1 = diffusion_vector,
  diffusion_vector.2 = diffusion_vector.2,
  summarisation = 'max',
  score_aggregation = 'min'
)

?tied_diffusion

?get_params

get_params(test_class, TRUE, TRUE)

test_class <- community_detection(
  test_class,
  min_seed_nodes = 0L,
  diffusion_threshold = .5,
  max_nodes = 300L
)

get_params(test_class, to_json = TRUE, pretty_json = TRUE)

x = get_results(test_class)

?get_results

test_class@params

devtools::load_all()

get_params(test_class, TRUE, TRUE)

?jsonlite::toJSON

jsonlite::toJSON(test_class@params) %>% jsonlite::prettify()

# RBH graph ----

protein_coding_genes <- data.table::fread("~/Desktop/protein_coding_genes.csv")

universe <- protein_coding_genes$id[1:500]

sets_per_origin <- 100

gene_sets_no <- sets_per_origin * length(LETTERS)

seed <- 123
set.seed(seed)

gene_sets <- purrr::map(1:gene_sets_no, ~ {
  set.seed(seed + .x + 1)
  size <- sample(20:100, 1)
  sample(universe, size, replace = FALSE)
})

names(gene_sets) <- purrr::map_chr(1:gene_sets_no, ~ {
  set.seed(seed + .x + 1)
  paste(sample(LETTERS, 5), collapse = "")
})

gene_sets_df <- purrr::imap(gene_sets, ~{
  data.table::data.table(
    genes = .x,
    name = .y
  )
})

origins <- rep(LETTERS, each = sets_per_origin)

gene_sets_df <- purrr::map2(gene_sets_df, origins, ~{
  .x[, origin := .y]
}) %>% rbindlist


module_df = gene_sets_df

devtools::document()

rbh_class = rbh_graph(
  gene_sets_df,
  dataset_col = 'origin',
  module_col = 'name',
  value_col = 'genes'
)

rbh_class



print(rbh_class)

rbh_class = generate_rbh_graph(rbh_class, minimum_similarity = .2, overlap_coefficient = F)

get_params(rbh_class, TRUE, TRUE)

list_of_list <- split(module_df %>% dplyr::select(!!module_col, !!value_col),
                      module_df[, ..dataset_col]) %>%
  purrr::map(., ~ {
    df <- .
    split(unlist(df[, ..value_col]),
          unlist(df[, ..module_col]))
  })


tictoc::tic()
rbh_results_v1 = rs_rbh_sets(list_of_list, FALSE, min_similarity = 0.2, FALSE)
tictoc::toc()

names(rbh_results_v1)

do.call(cbind, rbh_results_v1[c("origin_modules", "target_modules", "similarity")])

origin_vector <- unlist(purrr::map2(rbh_results_v1$origin, rbh_results_v1$comparisons, ~ {
  rep(.x, each = .y)
}))

target_vector <- unlist(purrr::map2(rbh_results_v1$target, rbh_results_v1$comparisons, ~ {
  rep(.x, each = .y)
}))

rbh_results_v1$origin_modules[rbh_results_v1$origin_modules == "NA"] <- NA

rbh_results_dt = data.table::data.table(
  origin = origin_vector,
  target = target_vector,
  origin_modules = rbh_results_v1$origin_modules,
  target_modules = rbh_results_v1$target_modules,
  similiarity = rbh_results_v1$similarity
) %>%
  .[!is.na(origin_modules) & similiarity >= 0.2] %>%
  .[, `:=`(
    combined_origin = paste(origin, origin_modules, sep = "_"),
    combined_target = paste(target, target_modules, sep = "_")
  )]

edge_df <- rbh_results_dt[, c("combined_origin", "combined_target", "similiarity")] %>%
  data.table::setnames(
    .,
    old = c("combined_origin", "combined_target", "similiarity"),
    new = c("from", "to", "weight")
  )

graph <- igraph::graph_from_data_frame(edge_df, directed = FALSE)


rbh_results_dt
