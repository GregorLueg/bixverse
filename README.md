# *bixverse package*

![r_package](https://img.shields.io/badge/R_package-0.0.1.4-orange) 

</br>

<img src="/misc/pics/bixverse_logo.png" width="128" height="128" alt="bixverse logo">

</br>

## *Description* 

This package contains various bioinformatics and computational biology workflows
that are being routinely used, ranging from gene set enrichment analyses, to 
network-based approaches for module detection in bulk RNAseq. The package provides
useful, bare bone versions of most bioinformatics functionalities and leverages Rust
to make any computational bottlenecks go *brrrrrrr* (i.e., fast).

## *Release notes*

We have now officially released version **0.0.1.4**. Test coverage has further increased,
and more and more functionality has been added. New features:

- CoReMo-based clustering, based on [Srivastava et al.](https://www.nature.com/articles/s41467-018-06008-4).
- TOM calculations, see [Yip and Horvath](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-22)
- Multi-level version of gsea also implemented for the gene ontology elimination
method.

## *Installation*

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
## *Docs*

- [Roadmap](/docs/roadmap.md)
- [Change log](/docs/change_log.md)
- [Why Rust](/docs/why_rust.md)
- [Code Style](/docs/code_style.md) (If you want to contribute).

## *Aim*

<img src="/misc/pics/but_why.png" width="418" height="218" alt="but why">

Good question, why this package? Basically, there are three reasons for this:
1. Initially, the package was born out of the desire to harmonise in a single 
package lots and lots of different packages. (One package to rule them all.) 
And while doing so, slim the functionality down to the bare bone essentials and
remove barely used options.
2. After having reimplemented for the 3rd time some version of a 
hypergeometric test for some project, one is just tired of it. The same applies 
to other established methods in the field that need some internal reimplementation
(with maybe some tweaks here and there) for reasons... ? Better to do it once
properly and make it public.
3. Rust makes everything in R so much faster, that we just wished to share the 
joy.

*Last update to the read-me: 16.06.2025*
