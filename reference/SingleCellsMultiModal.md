# bixverse SingleCells (multi modal) class

This is the `bixverse`-based SingleCells class for multiple modalities.
Under the hood it uses a DuckDB for obs and vars storing, and a
Rust-based binarised file format to store the raw and normalised counts
for single cell RNAseq. In both cases, the idea is not to hold any data
that is not needed at a given point of time in memory, but leverage
speedy on-disk computations and streaming engines powered by Rust and
DuckDB to run the analysis.

## Usage

``` r
SingleCellsMultiModal(dir_data)
```

## Arguments

- dir_data:

  String. This is the directory in which the experimental files will be

## Value

Returns the `SingleCellsMultiModal` class for further operations.

## Properties

- db_connection:

  This contains an R6 class with DuckDB pointers and wrappers to
  interact with the table-like data for this experiment.

- count_connection:

  This contains an R6-like environment that points to Rust functions
  that can work on the RNAseq counts more specifically.

- adt_counts:

  ...

- peak_connection:

  Future feature: to store ATAC Seq counts in the future.

- dir_data:

  Path to the directory in which the data will be saved on disk.

- sc_cache:

  Class with cached data. Contains less memory-heavy objects such as
  embeddings, kNN information or sNN graphs for the single cell RNAseq.

- adt_cache:

  Class with cached data. Contains less memory-heavy objects such as
  embeddings, kNN information or sNN graphs for the single cell
  Antibody-Derived Tags.

- atac_cache:

  Class with cached data. Contains less memory-heavy objects such as
  embeddings, kNN information or sNN graphs for the single cell
  chromatin accessability.

- sc_map:

  Class containing various mapping information such as HVG indices,
  cells to keep, etc.

- other_data:

  List that contains additional data and results, such as for example
  the WNN graph.

- dims:

  Dimensions of the original data.
