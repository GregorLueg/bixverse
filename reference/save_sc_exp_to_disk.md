# Save memory-bound data to disk

Helper function that stores the memory-bound data to disk for
checkpointing or when you close the session for quick recovery of prior
work. You have the option to save as `".rds"` or `".qs2"` (you need to
have the package `"qs2"` installed for this option!).

## Usage

``` r
save_sc_exp_to_disk(object, type = c("qs2", "rds"))
```

## Arguments

- object:

  `single_cell_exp` class.

- type:

  String. One of `c("qs2", "rds")`. Defines which binary format to use.
  Will default to `"qs2"` for speed.

## Value

The object with added information on the data on disk.
