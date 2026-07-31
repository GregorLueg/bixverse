# Constructor for a single pipeline step

Internal constructor used by the public `step_*()` functions. Wraps the
generic to dispatch on the `SingleCells`/`SingleCellsSubset` object
together with the captured arguments and a short human-readable name.

## Usage

``` r
new_sc_step(
  name,
  fn,
  args,
  accepts = c("SingleCells", "SingleCellsSubset"),
  returns = "input"
)
```

## Arguments

- name:

  String. Short identifier shown when printing the pipeline.

- fn:

  Function. Typically an S7 generic such as
  [`find_hvg_sc()`](https://gregorlueg.github.io/bixverse/reference/find_hvg_sc.md)
  that takes the object as the first argument.

- args:

  Named list. Arguments passed to `fn` at apply time (object is
  prepended automatically).

- accepts:

  Character vector. Class names the underlying generic has methods for.
  Checked by
  [`validate_pipeline()`](https://gregorlueg.github.io/bixverse/reference/validate_pipeline.md)
  before anything runs.

- returns:

  String. Class name the step returns, or `"input"` if it hands back
  whatever it was given.

## Value

An `ScStep` object.
