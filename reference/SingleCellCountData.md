# Single cell count data handler

**\[experimental\]** A class for handling single cell count data stored
on disk in two complementary binary representations: a CSR-like layout
(`f_path_cells`) for fast cell-wise access and a CSC-like layout
(`f_path_genes`) for fast gene-wise access. Both raw counts and
log-normalised counts are stored side by side. Provides methods for
ingesting data from R, `h5ad`, and `mtx` sources (including multi-file
workflows), converting between layouts, retrieving slices of the matrix,
and merging existing binary objects.

## Arguments

- f_path_cells:

  (`character`)  
  Path to the `.bin` file for the cell-based (CSR-like) representation.

- f_path_genes:

  (`character`)  
  Path to the `.bin` file for the gene-based (CSC-like) representation.

- n_cells:

  (`integer`)  
  Number of cells represented in the data.

- n_genes:

  (`integer`)  
  Number of genes represented in the data.

## Value

A new instance of the `SingleCellCountData` class.

## Methods

### Method `new`

Create a new instance of the class

#### Arguments

- `f_path_cells`:

  (`character`)  
  Path to the `.bin` file for the cell-based representation.

- `f_path_genes`:

  (`character`)  
  Path to the `.bin` file for the gene-based representation.

#### return

A new `SingleCellCountData` instance with `n_cells` and `n_genes`
initialised to zero.

### Method `get_shape`

Get the shape of the matrix

#### return

An integer vector `c(n_cells, n_genes)`.

### Method `set_from_file`

Populate `n_cells` and `n_genes` from the cells binary file

#### description

Reads the header of the file at `f_path_cells` and updates the `n_cells`
and `n_genes` fields accordingly. Useful when reconnecting to an
existing object on disk.

#### return

Invisible `NULL`.

### Method `r_data_to_file`

Write a CSR matrix from R to the cells binary file

#### Arguments

- `r_data`:

  (`list`)  
  A list convertible into `CompressedSparseData2`. Must contain the
  elements `"indptr"`, `"indices"`, `"data"`, `"nrow"`, `"ncol"` and
  `"format"`.

- `qc_params`:

  (`list`)  
  Quality control parameters parseable into `MinCellQuality`.

- `verbose`:

  (`logical`)  
  Controls verbosity of the function.

#### description

Ingest a sparse matrix passed in from R, apply per-cell QC, and write
the result to `f_path_cells`.

#### return

A list with `cell_indices`, `gene_indices`, `lib_size` and `nnz`.

### Method `h5ad_to_file`

Write an h5ad file to the cells binary file

#### Arguments

- `cs_type`:

  (`character`)  
  Storage layout of the h5ad data. One of `"CSC"` or `"CSR"`.

- `h5_path`:

  (`character`)  
  Path to the h5ad file.

- `no_cells`:

  (`integer`)  
  Number of cells in the h5ad file.

- `no_genes`:

  (`integer`)  
  Number of genes in the h5ad file.

- `qc_params`:

  (`list`)  
  Quality control parameters parseable into `MinCellQuality`.

- `slot`:

  (`character`)  
  Where to find the raw counts. One of `"X"` or `"raw.X"` (for CellXGene
  data).

- `verbose`:

  (`logical`)  
  Controls verbosity of the function.

#### return

A list with `cell_indices`, `gene_indices`, `lib_size` and `nnz`.

### Method `norm_h5ad_to_file`

Write an h5ad file with normalised counts to the cells binary file

#### Arguments

- `cs_type`:

  (`character`)  
  Storage layout of the h5 data. One of `"CSC"` or `"CSR"`.

- `h5_path`:

  (`character`)  
  Path to the h5 file.

- `no_cells`:

  (`integer`)  
  Number of cells in the h5 file.

- `no_genes`:

  (`integer`)  
  Number of genes in the h5 file.

- `obs_lib_size_col`:

  (`character`)  
  Name of the `obs` column containing total counts per cell (e.g.
  `"nCount_RNA"`).

- `target_size`:

  (`numeric`)  
  Target size used in the original normalisation (e.g. `1e4`).

- `qc_params`:

  (`list`)  
  Quality control parameters parseable into `MinCellQuality`.

- `verbose`:

  (`logical`)  
  Controls verbosity of the function.

#### description

For data sets where only normalised counts are available in `X`. Reads
library sizes from a specified `obs` column to reconstruct raw counts
before writing.

#### return

A list with `cell_indices`, `gene_indices`, `lib_size` and `nnz`.

### Method `h5ad_to_file_streaming`

Write an h5ad file to disk using streaming

#### Arguments

- `cs_type`:

  (`character`)  
  Storage layout of the h5 data. One of `"CSC"` or `"CSR"`.

- `h5_path`:

  (`character`)  
  Path to the h5 file.

- `no_cells`:

  (`integer`)  
  Number of cells in the h5 file.

- `no_genes`:

  (`integer`)  
  Number of genes in the h5 file.

- `qc_params`:

  (`list`)  
  Quality control parameters parseable into `MinCellQuality`.

- `slot`:

  (`character`)  
  Where to find the raw counts. One of `"X"` or `"raw.X"` (for CellXGene
  data).

- `verbose`:

  (`logical`)  
  Controls verbosity of the function.

#### description

Slower but lighter on memory than `h5ad_to_file`; streams the input
where possible.

#### return

A list with `cell_indices`, `gene_indices`, `lib_size` and `nnz`.

### Method `multi_h5ad_to_file`

Load multiple h5ad files into a single binary

#### Arguments

- `file_tasks`:

  (`list`)  
  A list of lists, each produced by the R prescan function. Each inner
  list must contain `exp_id`, `h5_path`, `cs_type`, `no_cells`,
  `no_genes` and `gene_local_to_universe`.

- `universe_size`:

  (`integer`)  
  Total number of genes in the universe.

- `qc_params`:

  (`list`)  
  Quality control parameters (`min_unique_genes`, `min_lib_size`,
  `min_cells`, `target_size`).

- `verbose`:

  (`logical`)  
  Controls verbosity.

#### return

A list with `global_gene_indices`, `total_cells`, `total_genes` and
`per_file` (a list of lists with `exp_id`, `cell_indices`, `lib_size`,
`nnz`).

### Method `mtx_to_file`

Write an mtx file to the cells binary file

#### Arguments

- `mtx_path`:

  (`character`)  
  Path to the mtx file.

- `qc_params`:

  (`list`)  
  Quality control parameters parseable into `MinCellQuality`.

- `cells_as_rows`:

  (`logical`)  
  `TRUE` if cells are rows in the mtx file, `FALSE` if cells are
  columns.

- `verbose`:

  (`logical`)  
  Controls verbosity of the function.

#### return

A list with `cell_indices`, `gene_indices`, `lib_size` and `nnz`.

### Method `mtx_to_file_streaming`

Write an mtx file to the cells binary file using streaming

#### Arguments

- `mtx_path`:

  (`character`)  
  Path to the mtx file.

- `qc_params`:

  (`list`)  
  Quality control parameters parseable into `MinCellQuality`.

- `cells_as_rows`:

  (`logical`)  
  `TRUE` if cells are rows in the mtx file, `FALSE` if cells are
  columns.

- `verbose`:

  (`logical`)  
  Controls verbosity of the function.

#### return

A list with `cell_indices`, `gene_indices`, `lib_size` and `nnz`.

### Method `multi_mtx_to_file`

Load multiple mtx files into a single binary

#### Arguments

- `file_tasks`:

  (`list`)  
  A list of lists, each containing `exp_id`, `mtx_path`, `cells_as_rows`
  and `gene_local_to_universe` (integer vector, `NA` for unmapped
  genes).

- `universe_size`:

  (`integer`)  
  Number of genes in the intersection universe.

- `qc_params`:

  (`list`)  
  Quality control parameters parseable into `MinCellQuality`.

- `verbose`:

  (`logical`)  
  Controls verbosity.

#### return

A list with `global_gene_indices`, `total_cells`, `total_genes` and
`per_file` (a list of lists with `exp_id`, `cell_indices`, `lib_size`,
`nnz`).

### Method `tenx_h5_to_file_streaming`

Write a 10x CellRanger h5 file to the cells binary file using streaming

#### Arguments

- `h5_path`:

  (`character`)  
  Path to the 10x h5 file.

- `version`:

  (`character`)  
  One of `"auto"`, `"v2"` or `"v3"`. `"auto"` detects the layout from
  the file.

- `no_cells`:

  (`integer`)  
  Number of cells (columns) in the file.

- `no_genes`:

  (`integer`)  
  Number of features (rows), including all modalities.

- `qc_params`:

  (`list`)  
  Quality control parameters parseable into `MinCellQuality`.

- `feature_type`:

  (`character` or `NULL`)  
  Target modality for v3. Defaults to `"Gene Expression"` when `NULL`.

- `verbose`:

  (`logical`)  
  Controls verbosity of the function.

#### description

Ingests the gene-expression modality from a CellRanger v2/v3 h5 file.
Other modalities (e.g. Antibody Capture) are filtered out via
`feature_type`.

#### return

A list with `cell_indices`, `gene_indices`, `lib_size` and `nnz`.

### Method `multi_tenx_h5_to_file`

Load multiple 10x CellRanger h5 files into a single binary

#### Arguments

- `file_tasks`:

  (`list`)  
  A list of lists, each produced by the R prescan function. Each inner
  list must contain `exp_id`, `h5_path`, `version` (`"v2"` or `"v3"`),
  `no_cells`, `no_genes`, `gene_local_to_universe` (integer vector, `NA`
  for unmapped / non-gene features) and `feature_type` (optional string,
  defaults to `"Gene Expression"`).

- `universe_size`:

  (`integer`)  
  Total number of genes in the universe.

- `qc_params`:

  (`list`)  
  Quality control parameters (`min_unique_genes`, `min_lib_size`,
  `min_cells`, `target_size`).

- `verbose`:

  (`logical`)  
  Controls verbosity.

#### return

A list with `global_gene_indices`, `total_cells`, `total_genes` and
`per_file` (a list of lists with `exp_id`, `cell_indices`, `lib_size`,
`nnz`).

### Method `return_full_mat`

Return the full count matrix

#### Arguments

- `assay`:

  (`character`)  
  One of `"raw"` or `"norm"`. Selects whether raw counts or
  log-normalised counts are returned.

- `cell_based`:

  (`logical`)  
  If `TRUE`, the data is returned in CSR layout (cells as rows). If
  `FALSE`, the data is returned in CSC layout (genes as columns).

- `verbose`:

  (`logical`)  
  Controls verbosity of the function.

#### return

A list with `indptr`, `indices`, `data`, `no_cells` and `no_genes`,
parseable into a sparse matrix in R.

### Method `get_cells_by_indices`

Return cells by index positions

#### Arguments

- `indices`:

  (`integer`)  
  The cell indices to return (1-indexed).

- `assay`:

  (`character`)  
  One of `"raw"` or `"norm"`.

#### description

Leverages the CSR-stored data for fast cell retrieval.

#### return

A list with `indptr`, `indices`, `data`, `no_cells` and `no_genes`,
parseable into a CSR matrix in R.

### Method `generate_gene_based_data`

Generate gene-based data from the cells binary file

#### Arguments

- `verbose`:

  (`logical`)  
  Controls verbosity of the function.

#### description

Reads the `.bin` file at `f_path_cells` and writes a gene-friendly (CSC)
representation to `f_path_genes`. The conversion happens fully in memory
and may cause memory pressure on large data sets; see
`generate_gene_based_data_streaming` or
`generate_gene_based_data_memory_bounded` for lighter alternatives.

#### return

Invisible `NULL`.

### Method `generate_gene_based_data_streaming`

Generate gene-based data with streaming

#### Arguments

- `batch_size`:

  (`integer`)  
  Number of cells processed per batch. Larger values increase memory
  pressure but reduce overhead.

- `verbose`:

  (`logical`)  
  Controls verbosity of the function.

#### description

Builds the CSC representation directly without creating intermediate CSR
structures. Suitable for very large data sets where the all-in- memory
path is too costly.

#### return

Invisible `NULL`.

### Method `generate_gene_based_data_memory_bounded`

Generate gene-based data with memory-bounded accumulation

#### Arguments

- `max_genes_in_memory`:

  (`integer`)  
  Maximum number of genes to accumulate at once (e.g. `2000`).

- `cell_batch_size`:

  (`integer`)  
  Number of cells to process at once (e.g. `100000`).

- `verbose`:

  (`logical`)  
  Controls verbosity.

#### description

Processes genes in phases to cap memory usage. Each phase:

1.  reads all cells (unavoidable for CSC conversion);

2.  only accumulates data for genes in the current phase;

3.  writes those genes to disk;

4.  clears memory and moves to the next phase.

#### return

Invisible `NULL`.

### Method `get_genes_by_indices`

Return genes by index positions

#### Arguments

- `indices`:

  (`integer`)  
  The gene indices to return (1-indexed).

- `assay`:

  (`character`)  
  One of `"raw"` or `"norm"`.

#### description

Leverages the CSC-stored data for fast gene retrieval.

#### return

A list with `indptr`, `indices`, `data`, `no_cells` and `no_genes`,
parseable into a CSC matrix in R.

### Method `get_nnz_genes`

Get the number of cells expressing each gene

#### Arguments

- `gene_indices`:

  (`integer` or `NULL`)  
  Optional 1-indexed gene indices. If `NULL`, results are returned for
  all genes.

#### return

An integer vector of NNZ counts for the requested genes.

### Method `merge_sc_files`

Merge multiple existing bin files into the cells binary file

#### Arguments

- `merge_tasks`:

  (`list`)  
  A list of lists. Each inner list must contain `exp_id`,
  `bin_cells_path`, `cells_to_keep` (0-indexed integer vector) and
  `gene_local_to_universe` (integer vector, `-1` for genes absent from
  the universe).

- `universe_size`:

  (`integer`)  
  Number of genes in the intersection universe.

- `renormalise`:

  (`logical`)  
  If `TRUE`, recompute `data_norm` against `target_size` using each
  cell's surviving raw counts. If `FALSE`, pass `data_norm` through
  untouched; the caller must guarantee all inputs were normalised
  against the same `target_size`.

- `target_size`:

  (`numeric`)  
  Target library size for renormalisation. Ignored when
  `renormalise = FALSE`.

- `verbose`:

  (`logical`)  
  Controls verbosity.

#### return

A list with `total_cells`, `total_genes` and `per_file` (a list of lists
with `exp_id`, `lib_size`, `nnz`).
