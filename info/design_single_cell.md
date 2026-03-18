# Technical design choices for the single cell implementation

This file is an attempt to explain the technical decisions behind the way the
singe cell has been implemented in this crate. This is an opinionated piece,
so take what is written here 

## Why do we need another single cell framework ... ?

<img src="../man/figures/standards.png" alt="standards">

Why the hell would you write a single cell framework from scratch and not use
one of the established libraries? We have in R [Seurat](https://satijalab.org/seurat/),
in the BioConductor universe [single cell experiments](https://bioconductor.org/books/release/OSCA/);
in Python [ScanPy](https://scanpy.readthedocs.io/en/stable/) and even now
GPU-accelerated [RAPIDS-singlecell](https://rapids-singlecell.readthedocs.io/en/latest/)
leveraging cuda and Nvidia GPUs. In Rust, we have [SingleRust](https://singlerust.com)
which is porting ScanPy over... What the hell is the point you might ask? Well,
the design philosophy in this package is quite different from all the solutions
above - if better/worse you shall discover, but it is different. But let's 
start with the problem setting...

### Chapter 1: Single cell is just not performant...

Some quotes from some recent papers: </br>

<i>We benchmarked rapids-singlecell on 1 million cells from the 10x Genomics mouse 
brain dataset, running a standard workflow: preprocessing, normalization, HVG 
selection, dimensionality reduction (PCA, UMAP, t-SNE), neighborhood graph 
construction, and Leiden clustering (Methods). <b>On a 32-core workstation, the 
pipeline took over 52 minutes (Methods)</b></i> - 
[from Dicks, et al., arXive, 2026](https://arxiv.org/abs/2603.02402)

<i>For instance, a typical single-cell data analysis pipeline5—encompassing data 
integration, clustering, visualization, and differential expression analysis—
requires about 16 h to process half a million cells on a standard desktop. 
When the cell number slightly increases to 600,000, <b>the above pipeline can crash 
due to memory exceeding, even on a professional computing platform with 512 GB 
RAM</b></i> - [from Li, et al., Nat Comm, 2025](https://www.nature.com/articles/s41467-025-56424-6)

</br>
Anyone who worked with more than 500k cells can attest to that it is painful.
VERY painful. The question then here is... why... ? Why can frameworks like
[DuckDB](https://duckdb.org) and [polars](https://pola.rs) deal with massive
data sets without needing fat memory machines? You can easily analyse 100's
of million of rows on a local laptop with these two. Why cannot we do this 
with single cell?

### Chapter 2: Do not keep data in memory...

Let's do some math... Let's say we have some single cell experiment with 1m
cells. We assume that in average something like 1000 genes are detected per
cell. What actually happens to our memory? And let's take the worst case here:
R.

```
# raw counts
x     (values, double in R by default):  1B × 8 bytes = 8 GB
                       (int32 if lucky): 1B × 4 bytes = 4 GB
i     (row indices, int32):              1B × 4 bytes = 4 GB
p     (col pointers, int32):             (1M+1) × 4   ~ 4 MB

Total: ~12 GB (double) / ~8 GB (int32)

# norm counts
For sure double, so 12 GB

# scaled counts with 3k HVG
2k HVGs:  1M × 2,000 × 8 bytes = 16 GB

# embeddings
PCA:              1M × 50 × 8 = 400 MB
Batch correction: 1M × 50 × 8 = 400 MB

Total: ~800 MB

# kNN graph
indices (int32):    20M × 4 = 80 MB
data (double):      20M × 8 = 160 MB
indptr (int32):      1M × 4 =   4 MB

Total: ~244 MB

# sNN graph
Same as kNN graph ~244 MB
```

How much are we looking at here? For sure 40 GBs of memory already occupied. 
Then you have some horrible tibbles with 1m of rows and just imagine a bunch
of string columns which can add in nasty situations for sure another 250 to 500
MB. Now R has a lovely copy and modify mechanic and suddenly everything becomes
unbearable slow, memory hungry and you end up with bloated objects that store 
way too much stuff in memory. The question is... Should we do any of this? Why 
do we keep the raw and norm counts in memory to start with? 
[BPCells](https://github.com/bnprks/BPCells) is going in the right direction
here, but can we do even better ... ?

### Chapter 3: Python and R are slow...

<img src="../man/figures/python_r_slow.jpg" alt="Python and R are slow">

There. I said it. They are slow and let us stop pretending otherwise. But meh,
NumPy! But meh, Rccp! Yeah great. You are using code/libraries in low level 
languages as an argument against these two languages being slow. You are still 
paying an performance and memory penalty just due to the way these languages are 
designed (which makes them great for quick iterative work!). You can make them 
faster, you will struggle to make them FAST. You do not have tight control over 
memory allocations, if you are not careful and you roundtrip between the 
kernels, things can get substantially slowed down and your memory pressure 
increases massively. 

## What is your proposal then ... ?

### Let's not keep data in memory...

Taking what is written above, we can now come up with a solution to the 
problems. Ideally, we want an interface from an interpreted, dynamically typed
language, so going all in into Rust is just too heavyhanded. But we can then
take certain design decisions:

1. We do not keep data in memory if not needed. Most analyses leveraging counts
   only need subsets of cells and genes. We can also avoid holding metadata in
   memory via DuckDB. Data will be streamed from disk when needed and Rust's
   compiler will free memory once it is not needed anymore. In times, this
   is even enforced by `drop()` to help. Raw counts? Nope. Normalised counts?
   Nope. Scaled counts? Nope. We just removed 40 GB of unnecessary stuff from
   memory. Data can be loaded into memory in ms if need be.
2. We do not even load all data into memory to start with. In 99.99% of the 
   cases, we know that we do not want low quality cells with tiny library sizes
   and lack of unique features. We also basically always to the usual library
   size to target size with log transformation. That means while reading in
   h5ad files or mtx files, we can already pre-scan these, only write 
3. We accept that single cell is noisy and do not even bother with `f32` or even
   `f64` precision. In terms of raw counts, most of the time `u16` is enough. 
   The ML field has shown us that one can get away with `f16`, 
   which reduces memory pressure substantially (they quantise even more 
   aggressively). We also do not bother doing any calculations in R (or Python). 
   Rust is superior for speed and optimisation, we can leverage aggressive 
   multi-threading, we avoid unnecessary copy semantics, we control
   memory tightly. R in this case just serves as an interface and telling
   the underlying Rust code what to do. More to that in the next section.

But if you do not have the data in memory, is this not going to be slow? Well
no. We store the data in structures like this. One for cells:

```{Rust}
/// CsrCellChunk
///
/// This structure is designed to store the data of a single cell in a
/// CSR-like format optimised for rapid access on disk.
///
/// Raw counts are stored as a `RawCounts` enum which transparently handles
/// both u16 and u32 storage. Gene indices are stored as u32 to support
/// data sets with more than 65,535 features.
#[derive(Debug)]
pub struct CsrCellChunk {
    /// Raw counts for this cell, stored as either u16 or u32 depending on
    /// the data range
    pub data_raw: RawCounts,
    /// Vector of the norm counts of this cell. A lossy compression for f16 is
    /// applied.
    pub data_norm: Vec<F16>,
    /// Total library size/UMI counts of the cell.
    pub library_size: usize,
    /// Index positions of the genes (u32 to support >65k features)
    pub indices: Vec<u32>,
    /// Original index in the data
    pub original_index: usize,
    /// Flag if the cell should be kept. (Not used at the moment.)
    pub to_keep: bool,
}
```

One for genes:

```{Rust}
/// CscGeneChunk
///
/// This structure is designed to store the data of a single gene in a
/// CSC-like format optimised for rapid access on disk.
///
/// Raw counts are stored as a `RawCounts` enum (u16 or u32). Cell indices
/// are stored as u32.
#[derive(Encode, Decode, Serialize, Deserialize, Debug)]
pub struct CscGeneChunk {
    /// Raw counts per cell/spot for this gene, stored as either u16 or u32
    pub data_raw: RawCounts,
    /// Normalised counts per cell/spot for this gene. Lossy compression to
    /// f16 is applied.
    pub data_norm: Vec<F16>,
    /// Average expression of this gene in the data.
    pub avg_exp: F16,
    /// Number of cells expressing this gene
    pub nnz: usize,
    /// Indices of the cells expressing this gene
    pub indices: Vec<u32>,
    /// Original index from the data
    pub original_index: usize,
    /// Flag to indicate if gene shall be kept. Not in use at the moment.
    pub to_keep: bool,
}
```

But you are duplicating the data! Yep, indeed. But this allows extremely 
efficient retrieval of data. You want to do DGEs between two groups? Only the
cells of these two groups will be loaded in via `CsrCellChunk`s. You want do
iterate over genes to identify the HVGs? That happens via the `CscGeneChunk`.
Another advantage is that in quite a few analyses, we can just load in subsets
of data and stream over the data. This allows to reduce the memory fingerprint
MASSIVELY during analyses. We still use [lz4](https://en.wikipedia.org/wiki/LZ4_(compression_algorithm))
compression to reduce the disk finger print, but the core philosophy is: </br></br>

<b>Disk space is cheap; memory ain't</b></br></br>

The approach allows to very quickly load in raw counts or normalised counts,
pending on what is needed and leverage the two-layer approach here (heavily
used by tileDB) to also not have to duplicate indices and indptr.

### Let's not round trip into interpreted languages

Continuing from step 3, we will also avoid round trips to R completely. This
avoids R quirks like copy-on-modify. One of the key algorithms in single cell
are (approximate) kNN searches, so, `bixverse` leverages a HIGHLY optimised
library for this, written from scratch, see [ann-search-rs](https://crates.io/crates/ann-search-rs).
The typical approximate nearest neighbour library used in R (Annoy) actually
writes the index to disk (which was Spotify's initial implementation)... Why
we are doing this in R? I don't know. There is also a full library for the
beloved (however not uncontroversial) 2D embeddings available that also 
leverages aggressively optimised Rust, see 
[manifoldsR](https://github.com/GregorLueg/manifoldsR). 
Million cell UMAP? Can be done in mere minutes... Any approach gets aggressively
optimised specifically for single cell in mind? SCENIC? Let's do multi-target
predictions with smart clustering of the genes together... The tree-based
regression models are anyway just used to calculate importance values, so, let
us aggressively remove everything that is not needed, leverage feature
quantisation, the fact that single cell is notoriously sparse and implement
highly optimised versions that run on local compute faster.

## There is no free launch