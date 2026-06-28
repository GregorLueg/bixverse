# bixverse SingleCells class

This is the `bixverse`-based SingleCells class. Under the hood it uses a
DuckDB for obs and vars storing, and a Rust-based binarised file format
to store the raw and normalised counts. In both cases, the idea is not
to hold any data that is not needed at a given point of time in memory,
but leverage speedy on-disk computations and streaming engines powered
by Rust and DuckDB to run the analysis. This version is specifically
designed for single cell RNAseq. If you want to use multi-modal data,
please refer to
[`SingleCellsMultiModal()`](https://gregorlueg.github.io/bixverse/reference/SingleCellsMultiModal.md) -
this class can store multiple layers of 'omics.

## Usage

``` r
SingleCells(dir_data)
```

## Arguments

- dir_data:

  String. This is the directory in which the experimental files will be

## Value

Returns the `SingleCells` class for further operations.

## Properties

- db_connection:

  This contains an R6 class with DuckDB pointers and wrappers to
  interact with the table-like data for this experiment.

- count_connection:

  This contains an R6-like environment that points to Rust functions
  that can work on the counts more specifically.

- dir_data:

  Path to the directory in which the data will be saved on disk.

- sc_cache:

  Class with cached data. Contains less memory-heavy objects such as
  embeddings, kNN information or sNN graphs.

- sc_map:

  Class containing various mapping information such as HVG indices,
  cells to keep, etc.

- dims:

  Dimensions of the original data.
