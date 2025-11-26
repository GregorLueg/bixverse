# fix hotspot ------------------------------------------------------------------

devtools::load_all()

dir_data <- path.expand("~/Downloads/single_cell_data/")

sc_object <- single_cell_exp(
  dir_data = dir_data
)

sc_object <- load_existing(sc_object)

sc_object <- set_cells_to_keep(
  x = sc_object,
  cells_to_keep = sc_object[[]][
    high_quality_cells == TRUE,
    cell_id
  ]
)

hvg_params <- params_sc_hvg(
  method = "vst"
)

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_params = hvg_params,
  hvg_no = 2000L,
  .verbose = TRUE
)

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = 30L,
  randomised_svd = TRUE,
  .verbose = TRUE
)

tictoc::tic()
umap_uwot <- uwot::umap(
  n_neighbors = 15,
  get_pca_factors(sc_object),
  metric = "cosine",
  n_epochs = 500L,
  verbose = TRUE
)
tictoc::toc()

?uwot::umap2

umap_uwot_dt <- as.data.table(umap_uwot) %>%
  `colnames<-`(c("umap1", "umap2")) %>%
  .[,
    cell_line := unlist(sc_object[["propagation_results"]], use.names = FALSE)
  ]

ggplot(data = umap_uwot_dt, mapping = aes(x = umap1, y = umap2)) +
  geom_point(mapping = aes(col = cell_line), size = 0.25)

# rextendr::document()
# rextendr::clean()

tictoc::tic()
umap_bixverse_sdg <- rs_umap(
  embd = get_pca_factors(sc_object),
  n_dim = 2L,
  ann_type = "annoy",
  optim = "sgd",
  k = 10L,
  42L,
  TRUE
)
tictoc::toc()

umap_uwot_bixverse_sgd <- as.data.table(umap_bixverse_sdg) %>%
  `colnames<-`(c("umap1", "umap2")) %>%
  .[,
    `:=`(
      cell_line = unlist(sc_object[["propagation_results"]], use.names = FALSE),
      condition = unlist(sc_object[["sample_combined"]], use.names = FALSE)
    )
  ]

ggplot(data = umap_uwot_bixverse_sgd, mapping = aes(x = umap1, y = umap2)) +
  geom_point(mapping = aes(col = cell_line), size = 0.25)

ggplot(data = umap_uwot_bixverse_sgd, mapping = aes(x = umap1, y = umap2)) +
  geom_point(mapping = aes(col = condition), size = 0.25)


tictoc::tic()
umap_bixverse_adam <- rs_umap(
  embd = get_pca_factors(sc_object),
  n_dim = 2L,
  ann_type = "annoy",
  optim = "adam",
  k = 15L,
  42L,
  TRUE
)
tictoc::toc()

umap_uwot_bixverse_adam <- as.data.table(umap_bixverse_adam) %>%
  `colnames<-`(c("umap1", "umap2")) %>%
  .[,
    `:=`(
      cell_line = unlist(sc_object[["propagation_results"]], use.names = FALSE),
      condition = unlist(sc_object[["sample_combined"]], use.names = FALSE)
    )
  ]

ggplot(data = umap_uwot_bixverse_adam, mapping = aes(x = umap1, y = umap2)) +
  geom_point(mapping = aes(col = cell_line), size = 0.25)

ggplot(data = umap_uwot_bixverse_adam, mapping = aes(x = umap1, y = umap2)) +
  geom_point(mapping = aes(col = condition), size = 0.25)
