# Seed the consensus k-means from the run seed

The Rust side reads a `consensus_seed` key for the k-means step.
[`params_nmf_consensus()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_consensus.md)
deliberately does not expose it: one call having a restart seed and a
separate clustering seed is a trap, where setting `seed` leaves the
k-means pinned at its default. So it is filled in here from whatever
`seed` the caller gave.

## Usage

``` r
.inject_consensus_seed(nmf_consensus_params, seed)
```

## Arguments

- nmf_consensus_params:

  List. Output of
  [`params_nmf_consensus()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_consensus.md).

- seed:

  Integer. The run seed.

## Value

The parameter list with `consensus_seed` set.
