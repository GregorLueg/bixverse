# bixverse: Rust-accelerated standard bioinformatics, computational biology workflows

Your favourite methods used in bioinformatics and computational biology
made VERY fast via Rust. Contains various gene set enrichment methods,
graph- based methods, co-expression module detection and much more.
Additionally, has now a single cell analysis suite, enabling million
cell scale analysis on local compute.

## Options

- `bixverse.cache_check`:

  Controls how strictly the single cell classes police stale cached
  artefacts (PCA, embeddings, kNN, sNN) against the current cell set.
  Unset, the default, is two tier: the getters warn, while
  [`assert_sc_state()`](https://gregorlueg.github.io/bixverse/reference/assert_sc_state.md)
  errors inside the functions that hand cached indices to Rust.
  `"error"` promotes the getter warnings to errors too, `"warn"` spells
  out the default getter behaviour, and `"none"` disables both tiers.
  See
  [`get_sc_cache_status()`](https://gregorlueg.github.io/bixverse/reference/get_sc_cache_status.md).

## See also

Useful links:

- <https://gregorlueg.github.io/bixverse>

## Author

**Maintainer**: Gregor Lueg <gregorlueg@me.com>

Authors:

- Gregor Lueg <gregorlueg@me.com>

Other contributors:

- Liesbeth Francois <francois.liesbeth@gmail.com> \[contributor\]

- Andrew Skelton <Ajskelton@proton.me> \[contributor\]

- Grant Neilson <grant.neilson5@gmail.com> \[contributor\]
