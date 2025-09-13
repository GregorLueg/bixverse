# single cell methods ----------------------------------------------------------

## i/o -------------------------------------------------------------------------

### h5ad -----------------------------------------------------------------------

load_h5ad <- S7::new_generic(
  name = "load_h5ad",
  dispatch_args = "object",
  fun = function(
    object,
    h5_path,
    min_genes = 200L,
    min_cells = 10L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

S7::method(load_h5ad, single_cell_exp) <- function(
  object,
  h5_path,
  min_genes = 200L,
  min_cells = 10L,
  .verbose = TRUE
) {}

#### helpers -------------------------------------------------------------------

# s3://cosyne-processed-data/cosyne-scrnaseq/ex053/BGI_Ultima/BGI/all-well/DGE_filtered/anndata.h5ad
