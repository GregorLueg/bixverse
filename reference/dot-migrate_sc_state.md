# Shape migration for state loaded off disk

[`save_sc_exp_to_disk()`](https://gregorlueg.github.io/bixverse/reference/save_sc_exp_to_disk.md)
writes the `ScMap` and `ScCache` lists verbatim, so an object saved
before a slot was added comes back missing it. Reading a missing element
gives `NULL` and nothing fails at load, which means the breakage
surfaces later and far from its cause.

Rebuilds the list from the current constructor and overlays whatever the
saved one carried. Keys the constructor no longer knows about are
dropped with a warning rather than silently kept, since they would never
be read again.

## Usage

``` r
.migrate_sc_state(x, template, label)
```

## Arguments

- x:

  The list loaded from disk.

- template:

  A freshly constructed list of the same class.

- label:

  Short label used in the warning.

## Value

The loaded state in the current shape.
