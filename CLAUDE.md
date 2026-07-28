# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## What this is

`bixverse` is an R package that wraps bioinformatics and computational
biology methods (gene set enrichment, graphs, co-expression modules,
ontologies, and a full single cell suite) around a Rust core exposed via
[rextendr](https://github.com/extendr/rextendr)/extendr. No official
Windows support (the h5 toolchain cross-compilation is the blocker).

## Commands

Build and install (compiles the Rust static lib, then regenerates
`R/extendr-wrappers.R`):

``` sh
R CMD INSTALL --preclean .
# or in R
devtools::install()
```

Tests use [tinytest](https://github.com/markvanderloo/tinytest), tests
live in `inst/tinytest/` and need the package installed:

``` r

tinytest::test_package("bixverse")                        # everything
tinytest::run_test_file("inst/tinytest/test_gsea.R")      # single file
```

Docs, formatting, Rust:

``` sh
air format .                                              # R formatter, 80 cols
cargo clippy --manifest-path=src/rust/Cargo.toml
```

``` r

devtools::document()   # roxygen2 8.0.0 -> man/, NAMESPACE
```

Pre-commit is configured (`air` formatting, large-file and R-artefact
guards).

`NOT_CRAN=1` in the environment enables non-CRAN build flags;
`tools/config.R` templates `src/Makevars{,.win}` from the `.in` files at
configure time. Builds are always release profile.

## Architecture

### Two-layer Rust split

The heavy numerical code does **not** live in this repo. It sits in the
standalone [`bixverse-rs`](https://crates.io/crates/bixverse-rs) crate,
declared in `src/rust/Cargo.toml` with features `single-cell` +
`multi-modal`. Check whether that dependency is a published version or a
local `path = ...` override before building: branches under active
development often point at a sibling `bixverse-rs` checkout, which must
exist at that path for the build to work. Everything under
`src/rust/src/` is a thin extendr binding layer, one module per domain
(`base`, `data`, `enrichment`, `graph`, `meta_cell`, `methods`,
`ontology`, `single_cell`), registered in `src/rust/src/lib.rs` via
`extendr_module!`.

Algorithmic work belongs in `bixverse-rs`; this repo gets the binding
plus the R wrapper.

`R/extendr-wrappers.R` is **generated** during compilation
(`cargo run --bin document`, see `src/rust/document.rs`). Never edit it
by hand. Exposed Rust functions are prefixed `rs_*`; each should have an
R wrapper that does the input validation.

### R layer conventions

Files are grouped by prefix:

- `classes_*.R` - class definitions (S7 mostly, R6 where encapsulation
  is needed)
- `methods_*.R` - S7 generics and their methods for those classes
- `functions_*.R` - plain procedural functions
- `param_wrappers.R` / `param_constructors.R` - `params_*()` functions
  returning validated parameter lists
- `checkmate_extensions.R` - `checkX()` / `assertX()` pairs built with
  [`checkmate::makeAssertionFunction()`](https://mllg.github.io/checkmate/reference/makeAssertion.html),
  mostly validating those parameter lists

Every exported function validates its inputs with checkmate. Parameter
bundles are passed as lists built by `params_*()` and validated by the
matching `assert*Params()`.

`BixverseBaseClass` (`R/base_class.R`) is the S7 root providing
[`get_params()`](https://gregorlueg.github.io/bixverse/reference/get_params.md)
/
[`get_results()`](https://gregorlueg.github.io/bixverse/reference/get_results.md);
most analysis classes inherit from it. `R/classes_union.R` defines S7
unions (e.g. `ScOrMc = SingleCells | MetaCells`) so methods can be
shared. Small result objects use plain S3; their `print` methods must be
registered manually in `R/zzz.R` (S7/S3 interaction quirk).

### Single cell design

`SingleCells` (`R/classes_single_cell_.R`) holds no counts in memory. It
is a handle over a directory containing:

- `counts_cells.bin` / `counts_genes.bin` - Rust binary count format,
  accessed through the `SingleCellCountData` extendr object (an
  env-based class, not R6)
- `sc_duckdb.db` - obs/var tables via the `SingleCellDuckDB` R6 class
  (`R/classes_single_cell_db.R`)

Two lightweight S3 side-cars ride along in the S7 object: `ScMap`
(gene/cell name to index maps, cells-to-keep, HVG indices) and `ScCache`
(PCA factors, embeddings, kNN, sNN graph). `SingleCellsMultiModal`
extends this with additional assay layers (currently ADT).

Gotcha: `ScMap` stores **0-based** indices for Rust. Setters like
[`set_cells_to_keep()`](https://gregorlueg.github.io/bixverse/reference/set_cells_to_keep.md)
subtract 1 on the way in; anything reading those slots and handing them
back to R must add 1.

## Code style

`info/code_style.md` is the authoritative document. The essentials:

- KISS. Avoid deep inheritance and speculative abstraction.
- Document everything with roxygen2; assert every input with checkmate.
- `data.table` over tibble/data.frame. `purrr::map*` over the apply
  family.
- Anything beyond simple reshaping/aggregation goes to Rust.
- Avoid new R dependencies; Rust dependencies are fine.
- 80 column limit, enforced by `air`.
- No emojis anywhere.

Vignettes are quarto (`.qmd` in `vignettes/`), site config in
`_pkgdown.yml`.

Version bumps in `DESCRIPTION` on `main` auto-tag and cut a GitHub
release; keep `src/rust/Cargo.toml` version in sync.
