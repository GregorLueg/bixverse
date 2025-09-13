# single cell methods ----------------------------------------------------------

## i/o -------------------------------------------------------------------------

### h5ad -----------------------------------------------------------------------

load_h5ad <- S7::new_generic(
  name = "load_h5ad",
  dispatch_args = "object",
  fun = function(object, .verbose) {
    S7::S7_dispatch()
  }
)

S7::method(load_h5ad, single_cell_exp) <- function(object, .verbose) {}

#### helpers -------------------------------------------------------------------

# s3://cosyne-processed-data/cosyne-scrnaseq/ex053/BGI_Ultima/BGI/all-well/DGE_filtered/anndata.h5ad
