# Assesses CoReMo-based gene module stability

The function assesses the stability of the CoReMo modules, leveraging a
leave-one-out sample method. In each iteration, one of the samples is
left out and the clustering is redone. Subsequently, the Jaccard
similarity (as a surrogate for membership stability) is calculated for
each feature across the different resamplings and added to the
`final_modules` data.table.

## Usage

``` r
cor_module_coremo_stability(object, chunk_size = 15L, .verbose = TRUE)
```

## Arguments

- object:

  The class, see [`bulk_coexp()`](bulk_coexp.md).

- chunk_size:

  Integer. Chunk size in which to process the data. Defaults to `15L`,
  i.e., 15 samples are being processed in one go. You can use bigger
  values here, but be aware of memory pressure.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

The class with added data to the properties for subsequent usage.

## References

Srivastava, et al., Nat. Commun., 2018; Francois, Romagnolo, et al.,
Nat. Commun., 2024.
