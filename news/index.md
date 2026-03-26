# Changelog

## bixverse 0.3.0

The big one… First, shift of versioning + release of full single cell
analysis suite. See below for more details.

- Full refactor of the Rust code, so it lives in its own
  [crate](https://crates.io/crates/bixverse-rs)
- Single cell officially released:
  - Class with getters that leverages on-disk storage of counts and a
    streaming engine to enable memory-friendly analysis of millions of
    cells.
  - I/O helpers that can read mtx outputs, h5ad, h5ad with only
    normalised counts, multiple h5ad and conversion from Seurat
    implemented.
  - Fully functional suite of pre-processing functions: proportion of
    gene sets, proportion of Top N genes of overall counts, HVG
    detection, PCA (with sparse solving avoiding densification),
    multiple kNN backends for different use-cases, sNN generation and
    Leiden community detection.
  - Three different batch correction methods implemented: BBKNN, fastMNN
    and Harmony.
  - Various analysis methods for single cell implemented: AUCell,
    GeneModuleScoring, VISION, Hotspot, SCENIC and miloR
    implementations.
  - Metacell generation: hdWGCNA bootstrap approach, SuperCells and
    SEACells.
  - Large number of convenience functions.
- Sister package started with [plotting
  functions](https://github.com/GregorLueg/bixverse.plots).
- Package with very fast 2D embedding methods also ready as a sister
  package – check out
  [manifoldsR](https://github.com/GregorLueg/manifoldsR).

## bixverse 0.0.2.3

- RCisTarget implemented in Rust, see [Aebir et
  al.](https://www.nature.com/articles/nmeth.4463) for details.
- Bug fix: scaling in
  [`preprocess_bulk_coexp()`](https://gregorlueg.github.io/bixverse/reference/preprocess_bulk_coexp.md)
  caused a matrix transpose.

## bixverse 0.0.2.2

- Mitch multi contrast enrichment implemented in Rust.
- Rework and refactor of sub modules in the Rust code.
- Additional Rust functions for set similarity across a list of items in
  R.
- QC function for the bulk DGE class split into
  [`qc_bulk_dge()`](https://gregorlueg.github.io/bixverse/reference/qc_bulk_dge.md)
  and
  [`normalise_bulk_dge()`](https://gregorlueg.github.io/bixverse/reference/normalise_bulk_dge.md).
- Bulk DGE function split into
  [`calculate_dge_limma()`](https://gregorlueg.github.io/bixverse/reference/calculate_dge_limma.md)
  and
  [`calculate_dge_hedges()`](https://gregorlueg.github.io/bixverse/reference/calculate_dge_hedges.md).
- Splitting by correlation sign and Eigengene calculation added for
  CoReMo.
- SNF algorithm added from [Wang et
  al.](https://www.nature.com/articles/nmeth.2810).
- Spectral clustering methods implemented in Rust.
- Breaking change:
  [`rs_set_similarity_list()`](https://gregorlueg.github.io/bixverse/reference/rs_set_similarity_list.md)
  renamed to
  [`rs_set_similarity_list2()`](https://gregorlueg.github.io/bixverse/reference/rs_set_similarity_list2.md).
  The former now takes a single list and calculates set similarity
  across all permutations of elements.
- Breaking change:
  [`preprocess_bulk_dge()`](https://gregorlueg.github.io/bixverse/reference/preprocess_bulk_dge.md)
  and
  [`calculate_all_dges()`](https://gregorlueg.github.io/bixverse/reference/calculate_all_dges.md)
  have been removed.

## bixverse 0.0.2.1

- Constrained PageRank implementation added.
- Improved synthetic data for correlation-based gene module detection.
- Improved tests for ICA, CoReMo and graph-based correlation modules.
- Various distance calculations accelerated in Rust.
- Rust-based personalised PageRank implementations now support weighted
  graphs.
- Bug fixes for the graph-based correlation methods.
- Fixed export bugs for some plotting functions.
- Breaking change: interface to Rust-based personalised PageRank
  calculations may break due to the addition of a `weights` parameter.

## bixverse 0.0.2.0

- Rework of the ontology class with Wang similarity added as an
  additional measure.
- Rework of the genetic community detection class and diffusion methods,
  including permutation-based testing for diffusions.
- Rework of the ICA code with a testing suite added for the ICA class.
- Reciprocal best hit method based on correlations added; an additional
  parameter is now required specifying whether to use set similarity or
  correlation-based similarity.
- Vignettes added for various methods.
- GSVA and ssGSEA implemented in Rust.
- Rust code reworked in various places for improved performance.
- Mutual information calculations implemented in Rust.
- Improvements to the synthetic data.
- Bug fixes for the graph-based correlation methods.
- Breaking change: ontology class functions and methods have been
  modified heavily and may be breaking.
- Breaking change: community detection method has been updated; old code
  may not work.
- Breaking change: ICA detection functions and methods have been updated
  and may break existing code.

## bixverse 0.0.1.4

- TOM implemented in Rust for correlation-based module detection.
- CoReMo-based gene module detection added.
- Breaking change: `cor_module_check_res()` renamed to
  [`cor_module_graph_check_res()`](https://gregorlueg.github.io/bixverse/reference/cor_module_graph_check_res.md)
  and `cor_module_final_modules()` renamed to
  [`cor_module_graph_final_modules()`](https://gregorlueg.github.io/bixverse/reference/cor_module_graph_final_modules.md).

## bixverse 0.0.1.3

- fgsea multi level method implemented; `n_more_extreme` is now returned
  by GSEA functions.
- RBF function implementations added for R matrices.

## bixverse 0.0.1.2

- fgsea simple implemented for gene ontology with elimination method.
- Improved gene ontology elimination methods to reduce unnecessary
  copying.

## bixverse 0.0.1.1

- Semantic similarities for ontologies added.
- Speed improvements in various Rust functions via reduced memory
  copying and use of lifetimes.
- Package stability improved with tests powered by
  [tinytest](https://github.com/markvanderloo/tinytest).
- Set similarity Rust functions exposed.
- DGE class added for Limma Voom differential gene expression
  calculations.
- Wrapper functions added to load h5ad objects into R memory.
- Traditional GSEA and fgsea simple ported into Rust.
- Bug fix:
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  for iterating over resolutions in reciprocal best hit graph
  generation.
- Bug fix: hypergeometric calculations and RBH graph.
- Breaking change:
  [`community_detection()`](https://gregorlueg.github.io/bixverse/reference/community_detection.md)
  now uses
  [`params_community_detection()`](https://gregorlueg.github.io/bixverse/reference/params_community_detection.md)
  for parameters.

## bixverse 0.0.1.0

- Hypergeometric tests for gene set analysis for single or lists of
  target gene sets.
- Gene ontology aware hypergeometric test for GO analysis.
- Network diffusion methods based on personalised PageRank (single and
  tied diffusion).
- Community detection algorithms on top of diffusion scores.
- Reciprocal best hit graphs using set similarities between gene modules
  from different data sets or methods.
- Contrastive PCA implementation for gene module detection.
- Differential correlation-based methods for gene module detection using
  network-based graph community detection.
- Independent component analysis with stabilised versions.
- Helper functions for Hedge’s G effect size and OpenTargets score
  summarisation.
