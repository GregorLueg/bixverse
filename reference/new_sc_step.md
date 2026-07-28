# Constructor for a single pipeline step

Internal constructor used by the public `step_*()` functions. Wraps the
generic to dispatch on the `SingleCells`/`SingleCellsSubset` object
together with the captured arguments and a short human-readable name.

## Usage

``` r
new_sc_step(name, fn, args)
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

## Value

An `ScStep` object.
