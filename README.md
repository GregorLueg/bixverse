# bixverse package

![r_package](https://img.shields.io/badge/R_package-0.3.2-orange) 
[![CI](https://github.com/GregorLueg/bixverse/actions/workflows/R-cmd-check.yml/badge.svg)](https://github.com/GregorLueg/bixverse/actions/workflows/R-cmd-check.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![pkgdown](https://img.shields.io/badge/pkgdown-website-1b5e9f?logo=github)](https://gregorlueg.github.io/bixverse/)
[![extendr](https://img.shields.io/badge/extendr-^0.9.0-276DC2)](https://extendr.github.io/extendr/extendr_api/)

</br>

<img src="man/figures/bixverse_logo.png" width="128" height="128" alt="bixverse logo">

</br>

## Intro

### Description

This is an *opionated* package making various bioinformatics workflows in R (or
ported from Python) much faster via low-level implementations in Rust. The core 
idea is to take different methods, write implementations in a compiled,
memory-managed language with minimal kernel round trips and leverage R purely as 
an orchestraction layer.
Result? Blazingly fast performance with low memory usage, making large-scale
analyses feasable without any cloud compute. Over time more and more methods
will be added. The aim will be to come a `tidyverse` equivalent, but for
a lot of downstream methods post WGS processing (with a strong emphasis on 
transcriptomics, think single cell and RNASeq for now). There is a sister 
package for plotting functions being build in parallel, see 
[here](https://github.com/GregorLueg/bixverse.plots) (that one is in alpha
phase).

<img src="/man/figures/bixverse_single_cell.png" width="500" height="500" alt="atlas scale on small compute">

### Release notes

Major release with `"0.4.0". The single cell suite has been further improved. 
The key changes are:
- Performance improvements across several axis. A lot of the underlying Rust
  code was made even faster and more performant. The goal of analysing a 1m
  single cell data set on a MacBook Air with 24 GB memory has been achieved.
- Multi-modal support. There is now a `SingleCellsMultiModal` class that allows
  you to (for now only) also add ADT counts. Some readers have been added to
  load in the data. This also includes support for new methods in this space:
  * [DSB normalisation](https://www.nature.com/articles/s41467-022-29356-8) for ADT counts
  * [WNN generation](https://www.cell.com/cell/fulltext/S0092-8674(21)00583-3) 
    for generation of multi-modal graphs.
  * Large number of methods and updates to support multi-model analysis
- Massive improvements on the two sister packages ([bixverse.plots](https://github.com/GregorLueg/bixverse.plots) and 
  [bixverse.gpu](https://github.com/GregorLueg/bixverse.gpu)): enabling 
  GPU-accelerated methods (Harmony, PCA) and a large number of plotting helpers.
- Updates to various vignettes to reflect the changes with this release.

[Please checkout out the website of the package for details](https://gregorlueg.github.io/bixverse/) -> 
particularly the sections around single cell (design choices and vignettes.)

## Usage

### Installation

You will need Rust on your system to have the package working. An installation
guide is provided [here](https://www.rust-lang.org/tools/install). There is a
bunch of further help written [here](https://extendr.github.io/rextendr/index.html)
by the rextendr guys in terms of Rust set up. (bixverse uses rextendr to interface
with Rust.)

Steps for installation:

1. In the terminal, install [Rust](https://www.rust-lang.org/tools/install)

```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

2. In R, install [rextendr](https://extendr.github.io/rextendr/index.html):

```
install.packages("rextendr")
```

3. Finally install bixverse:

```
devtools::install_github("https://github.com/GregorLueg/bixverse")
```

### Windows support

If you are using Windows, I am sorry, the tool chain is just very, very 
painful... I really tried to make this work and maybe there are some hacks in 
terms of compiling everything to install the package, but it has proven...
challenging in the CI/CD. Hence, no official Windows support for now. It is
specifically the incorporation of h5 which proves non-trivial with 
cross-compiling that with Rust within the R umbrella.

### How to use the package.

The package website can be found [here](https://gregorlueg.github.io/bixverse/).
A good primer for [why Rust is here](https://gregorlueg.github.io/bixverse/articles/rust_functions.html) - a show case of how much faster Rust can make a lot of basic functions much faster. 
If you wish to integrate this into your package, please feel free. If
you wish to use the single cell part, it is really worth reading this 
[here](https://gregorlueg.github.io/bixverse/articles/design_single_cell.html) 
first... It will give you a good explainer on the design decisions, the choices
and trade-offs. The various vignettes will show you how to analyse data.

## Roadmap

### Single and spatial transcriptomics:

#### General

- [x] Multi h5 (10x output) i/o
- [ ] Save data to h5ad for easier interoperability with Python.
- [ ] Splitting a SingleCells into several sub directories for easier management
- [ ] Implementation of [Palantir](https://www.nature.com/articles/s41587-019-0068-4) and
  [Slingshot](https://pubmed.ncbi.nlm.nih.gov/29914354/) for trajectory
  analysis
- [ ] Port over [NicheNet](https://www.nature.com/articles/s41592-019-0667-5).
- [ ] Easy interoperability that chunks of data can be read in for neural 
  network training in the corresponding deep learning frameworks.

#### Spatial 

Full support of spatial transcriptomics via a dedicated `SpatialSpots` class
leveraging the current Rust-based infrastructure for large scale analysis on
local compute.

### General methods

- [ ] Addition of an NMF method for dense and sparse matrices. After some 
  research likely a combination HALS + NNDSVD initialisation

### Further cross integration with other packages

Two sister packages are in active development: 
[GPU-accelerated methods](https://github.com/GregorLueg/bixverse.gpu) and a
dedicated [plotting package](https://github.com/GregorLueg/bixverse.plots). 
Further improved cross integration is on the to-do list.

### Documentation

While there are already quite a few vignettes, the amount of code in the package
is... quite substantial and there are methods hidden here and there that lack
any vignettes for now. On the roadmap to provide examples here, too.

## For developers

If you wish to contribute, please read the [Code Style](/info/code_style.md).
