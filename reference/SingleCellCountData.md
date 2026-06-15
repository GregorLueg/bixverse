# A class for handling single cell count data

### Params

## Usage

``` r
SingleCellCountData
```

## Details

- `f_path_cells` - Path to the .bin file for the cells.

- `f_path_genes` - Path to the .bin file for the genes.

- `n_cells` - No of cells represented in the data.

- `n_genes` - No of genes represented in the data.

## Methods

### Method `new`

Create new instance of the class

#### Params

- `f_path_cells` - Path to the .bin file for the cells.

- `f_path_genes` - Path to the .bin file for the genes.

#### Method `get_shape`

Get the shape

#### Returns

Vector with rows x cells

#### Method `r_data_to_file`

Write data from R CSR to disk

Helper function to write CSR matrices from R to disk.

#### Params

- `r_data` - A list that can be transformed into CompressedSparseData2.
  Needs to have following elements: `"indptr"`, `"indices"`, `"data"`,
  `"nrow"`, `"ncol"` and `"format"`.

- `qc_params` - List with the quality control parameters.

- `verbose` - Controls verbosity of the function.

#### Returns

A list with QC parameters.

#### Method `h5_to_file`

Save h5 to file

#### Params

- `cs_type` - How was the h5 data saved. CSC or CSR.

- `h5_path` - Path to the h5 file.

- `no_cells` - Number of cells in the h5 file.

- `no_genes` - Number of genes in the h5 file.

- `qc_params` - List with the quality control parameters.

- `verbose` - Controls verbosity of the function.

#### Returns

A list with qc parameters.

#### Method `norm_h5_to_file`

Save h5 with normalised counts to file

For datasets where only normalised counts are available in X. Reads
library sizes from a specified obs column to reconstruct raw counts.

#### Params

- `cs_type` - How was the h5 data saved. CSC or CSR.

- `h5_path` - Path to the h5 file.

- `no_cells` - Number of cells in the h5 file.

- `no_genes` - Number of genes in the h5 file.

- `obs_lib_size_col` - Name of the obs column containing total counts
  per cell (e.g. "nCount_RNA").

- `target_size` - Target size used in the original normalisation (e.g.
  1e4).

- `qc_params` - List with the quality control parameters.

- `verbose` - Controls verbosity of the function.

#### Returns

A list with qc parameters.

#### Method `h5_to_file_streaming`

Save h5 to file

Slower version that is less memory heavy and will make usage of
streaming where possible.

#### Params

- `cs_type` - How was the h5 data saved. CSC or CSR.

- `h5_path` - Path to the h5 file.

- `no_cells` - Number of cells in the h5 file.

- `no_genes` - Number of genes in the h5 file.

- `qc_params` - List with the quality control parameters.

- `verbose` - Controls verbosity of the function.

#### Returns

A list with qc parameters.

#### Method `multi_h5_to_file`

Load multiple h5ad files into a single binary

#### Params

- `file_tasks` - R list of lists, each produced by the R prescan
  function. Each inner list must contain: exp_id, h5_path, cs_type,
  no_cells, no_genes, gene_local_to_universe.

- `universe_size` - Total number of genes in the universe.

- `qc_params` - List with QC parameters (min_unique_genes, min_lib_size,
  min_cells, target_size).

- `verbose` - Controls verbosity.

#### Returns

A list with: global_gene_indices, total_cells, total_genes, per_file
(list of lists with exp_id, cell_indices, lib_size, nnz).

#### Method `mtx_to_file`

Save mtx to file

#### Params

- `mtx_path` - Path to the mtx file.

- `qc_params` - List with the quality control parameters.

- `cells_as_rows` - Do the cells represent rows (= TRUE) or columns.

- `verbose` - Controls verbosity of the function.

#### Returns

A list with qc parameters.

#### Method `mtx_to_file_streaming`

Save mtx to file - streaming version

#### Params

- `mtx_path` - Path to the mtx file.

- `qc_params` - List with the quality control parameters.

- `cells_as_rows` - Do the cells represent rows (= TRUE) or columns.

- `verbose` - Controls verbosity of the function.

#### Returns

A list with qc parameters.

#### Method `return_full_mat`

Returns the full matrix

#### Params

- `assay` - String. Return the raw counts or log-normalised counts. One
  of `"raw"` or `"norm"`.

- `cell_based` - Boolean. Shall the data be returned in CSR or CSC.

- `verbose` - Boolean. Verbosity of the function.

#### Returns

An R list with all of the info that was stored in the .bin file

#### Method `get_cells_by_indices`

Return cells by index positions

Leverages the CSR-stored data for fast cell retrieval

#### Params

- `indices` - The cell indices which to return (1-indexed).

- `assay` - Shall the raw or norm counts be returned

#### Returns

A list that can be parsed into a CSR matrix in R

#### Method `generate_gene_based_data`

Transforms already written cell data also into the gene data

This function will read the .bin file at `self.f_path_cells` and
transform the data into the gene-based file format. This happens for the
full data set in memory and might cause memory pressure, pending the
size of the data.

#### Params

- `qc_params` - List with the quality control parameters.

- `verbose` - Controls verbosity of the function.

#### Method `generate_gene_based_data_streaming`

Generate gene-based data with streaming to reduce memory pressure

This approach builds CSC format directly without creating intermediate
CSR structures. Ideal for very large data sets.

#### Params

- `batch_size` - Size of the batch to process in one go. The larger, the
  more memory pressure will occur.

- `verbose` - Controls verbosity of the function.

#### Method `generate_gene_based_data_memory_bounded`

Generate gene-based data with memory-bounded accumulation

This processes genes in phases to limit memory usage. Each phase:

1.  Reads all cells (unavoidable for CSC conversion)

2.  Only accumulates data for genes in current phase

3.  Writes those genes to disk

4.  Clears memory and moves to next phase

#### Params

- `max_genes_in_memory` - Maximum genes to accumulate at once (e.g.,
  2000)

- `cell_batch_size` - How many cells to process at once (e.g., 100000)

- `verbose` - Controls verbosity

#### Method `get_genes_by_indices`

Return genes by index positions

Leverages the CSC-stored data for fast gene retrieval

#### Params

- `indices` - The gene indices which to return (1-indexed).

- `assay` - Shall the raw or norm counts be returned

#### Returns

A list that can be parsed into a CSC matrix in R

#### Method `set_from_file`

Set cell numbers and genes

#### Params

- `cell_no` - No of cells

- `gene_no` - No of genes

#### Method `get_nnz_genes`

Helper function to get the number of cells expressing a gene

#### Params

- `gene_indices` - Optional gene indices (1-indexed!). If None, returns
  all genes.

#### Returns

A vector of NNZ for the genes.
