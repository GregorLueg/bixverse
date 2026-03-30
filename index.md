# bixverse package

![r_package](https://img.shields.io/badge/R_package-0.3.0-orange)[![CI](https://github.com/GregorLueg/bixverse/actions/workflows/R-cmd-check.yml/badge.svg)](https://github.com/GregorLueg/bixverse/actions/workflows/R-cmd-check.yml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![pkgdown](https://img.shields.io/badge/pkgdown-website-1b5e9f?logo=github)](https://gregorlueg.github.io/bixverse/)

![bixverse logo](reference/figures/bixverse_logo.png)

## Intro

### Description

This is an *opionated* package making various bioinformatics workflows
in R (or ported from Python) much faster via low-level implementations
in Rust. The core idea is to take different methods, write
implementations in a compiled, memory-managed language with minimal
kernel round trips and leverage R purely as an orchestraction layer.
Result? Blazingly fast performance with low memory usage, making
large-scale analyses feasable without any cloud compute. Over time more
and more methods will be added. The aim will be to come a `tidyverse`
equivalent, but for a lot of downstream methods post WGS processing.
There is a sister package for plotting functions being build in
parallel, see [here](https://github.com/GregorLueg/bixverse.plots) (that
one is in alpha phase).

### Release notes

With the `0.3.0` a lot has happened. The lack of updates had a reason…
Some heavy cooking has been going on… The package now contains a full
release of the a single cell functionality suite that you can use to
analyse millions of cells locally and implements already a large number
of methods in this space into Rust and exposes them to R. [Please
checkout out the website of the
package.](https://gregorlueg.github.io/bixverse/), specifically the
sections arounds single cell.

## Usage

### Installation

You will need Rust on your system to have the package working. An
installation guide is provided
[here](https://www.rust-lang.org/tools/install). There is a bunch of
further help written
[here](https://extendr.github.io/rextendr/index.html) by the rextendr
guys in terms of Rust set up. (bixverse uses rextendr to interface with
Rust.)

Steps for installation:

1.  In the terminal, install
    [Rust](https://www.rust-lang.org/tools/install)

&nbsp;

    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

2.  In R, install
    [rextendr](https://extendr.github.io/rextendr/index.html):

&nbsp;

    install.packages("rextendr")

3.  Finally install bixverse:

&nbsp;

    devtools::install_github("https://github.com/GregorLueg/bixverse")

### Windows support

If you are using Windows, I am sorry, the tool chain is just very, very
painful… I really tried to make this work and maybe there are some hacks
in terms of compiling everything to install the package, but it has
proven… challenging in the CI/CD. Hence, no official Windows support for
now. It is specifically the incorporation of h5 which proves non-trivial
with cross-compiling that with Rust within the R umbrella.

### How to use the package.

The package website can be found
[here](https://gregorlueg.github.io/bixverse/). A good primer for [why
Rust is
here](https://gregorlueg.github.io/bixverse/articles/rust_functions.html) -
a show case of how much faster Rust can make a lot of basic functions
much faster. If you wish to integrate this into your package, please
feel free. If you wish to use the single cell part, it is really worth
reading this
[here](https://gregorlueg.github.io/bixverse/articles/design_single_cell.html)
first… It will give you a good explainer on the design decisions, the
choices and trade-offs. The various vignettes will show you how to
analyse data.

## Roadmap

### Single and spatial transcriptomics:

- For single cell the following stuff will be hopefully soon’ish
  implemented:
  - More multi-file read in support. At the moment, multiple h5ad files
    are possible, but not yet for other file formats.
  - h5 file i/o (provided by CellRanger).
  - Saving data to h5ad for easier interoperability with Python.
  - Something cool with [Zarr](https://zarr.dev) … ?
  - Methods on top of the meta cells: co-expression network detection
    etc.
  - Expansion of the [sister
    package](https://github.com/GregorLueg/bixverse.plots) to have
    plotting helpers in there for single cell.
  - Helpers to slice and dice the data easier and add new data - this
    will
  - Implementations of
    [Palantir](https://www.nature.com/articles/s41587-019-0068-4) and
    [Slingshot](https://pubmed.ncbi.nlm.nih.gov/29914354/).
  - Port over
    [NicheNet](https://www.nature.com/articles/s41592-019-0667-5)
  - Add more GPU-acceleration via [cubecl/WGPU backend](NA) for
    GPU-agnostic acceleration where appropriate, see another [sister
    package](https://github.com/GregorLueg/bixverse.gpu)
- Leverage the current infrastructure and add dedicated support and
  methods for spatial transcriptomics. There are some cool methods in
  that space that for sure could benefit from the speed that a compiled,
  memory-managed language offers. Especially when analysing more data
  sets.
- Add other interesting methods that I can find (and have a use-case
  for).

### Documentation

While there are already quite a few vignettes, the amount of code in the
package is… quite substantial and there are methods hidden here and
there that lack any vignettes for now. On the roadmap to provide
examples here, too.

## For developers

If you wish to contribute, please read the [Code
Style](https://gregorlueg.github.io/info/code_style.md).
