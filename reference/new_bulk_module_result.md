# Constructor for a bulk co-expression module result

Uniform container for the terminal output of every bulk co-expression
method. Stores the module membership (gene to module_id assignment), the
method-specific factor matrices (gene loadings, sample activities,
eigengenes, dictionaries), the fit parameters, and method-specific
diagnostics.

## Usage

``` r
new_bulk_module_result(
  modules,
  factors = list(),
  method,
  params = list(),
  diagnostics = list()
)
```

## Arguments

- modules:

  data.table. Module membership. Must contain `gene` and `module_id`
  columns. Method-specific extras (`sign`, `stability`, `loading`, ...)
  are allowed.

- factors:

  Named list of matrices. Method-agnostic keys shared across methods:
  `gene_loadings`, `sample_activity`, `module_eigengenes`, `dictionary`,
  `loadings`, `gene_to_eigengene_cor`. May be empty for methods that
  produce no factor matrices (e.g. Leiden).

- method:

  String. One of `c("correlation-based",`
  `"differential correlation-based", "ICA-based", "dgrdl-based",`
  `"nmf-based")`.

- params:

  List. The `params_xxx()` list that produced the fit.

- diagnostics:

  Named list. Method-specific diagnostics (stability data.table,
  resolution used, per-run losses, laplacians, ...).

## Value

An object of class `BulkModuleResult`.
