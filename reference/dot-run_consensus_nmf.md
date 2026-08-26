# Run a consensus NMF binding with an actionable error

The consensus step can only fail once the restarts are done: the density
filter may leave fewer than `k` components, or a cluster may end up
empty. Both are real signals that the factorisation is not reproducible
at this `k`, so neither is retried or degraded into a partial answer.
The Rust message says what happened; this appends what to do about it.

## Usage

``` r
.run_consensus_nmf(.rs_call, nmf_consensus_params, seed, ...)
```

## Arguments

- .rs_call:

  Function. The `rs_nmf_consensus_*` binding to call.

- nmf_consensus_params:

  List. Output of
  [`params_nmf_consensus()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_consensus.md).
  The consensus seed is injected here.

- seed:

  Integer. The run seed.

- ...:

  Remaining arguments forwarded to `.rs_call`.

## Value

The raw list the binding returned.
