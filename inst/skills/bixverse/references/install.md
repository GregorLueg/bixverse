# Installing bixverse and its ecosystem

## The package itself

bixverse compiles a Rust static library at install time, so a Rust toolchain has
to be on the system first. Three steps, in order.

1. Rust, in the terminal:

```sh
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

2. rextendr, in R:

```r
install.packages("rextendr")
```

3. bixverse:

```r
devtools::install_github("https://github.com/GregorLueg/bixverse")
```

The first install compiles the whole Rust dependency tree and takes a while.
That's normal, not a hang.

`SystemRequirements` is `Cargo` and `rustc`. R 4.2 or newer.

## Windows

Not supported. The blocker is incorporating h5 when cross-compiling Rust under
the R umbrella, and it does not work in CI. Hacks may exist but nothing is
supported, so don't burn time inventing one. Use WSL, macOS or Linux.

## Optional dependencies

Most of bixverse works out of the box. These `Suggests` are pulled in only by
specific methods, and a missing one gives an error naming the package:

| Package | Needed for |
|---|---|
| `bixverse.plots` | nearly all single cell plotting |
| `Seurat` | `load_seurat()` |
| `fgsea` | benchmarking against the reference fgsea implementation |
| `GSVA` | benchmarking against reference GSVA |
| `singscore`, `mitch` | reference implementations for the pathway activity methods |
| `msigdbr` | pulling MSigDB gene sets |
| `biomaRt` | gene identifier conversion |
| `ontologyIndex`, `ontologySimilarity` | ontology cross-checks |
| `qs2` | `save_sc_exp_to_disk(type = "qs2")`, the faster serialisation path |
| `quarto` | building vignettes |

Note that Scrublet and scDblFinder are **reimplemented in Rust** inside
bixverse. `scrublet_sc()` and `scdblfinder_sc()` do not need the original
Python or Bioconductor packages.

## The sister packages

bixverse is deliberately narrow. Several things a user will ask for live next
door.

| Package | What it is | Reach for it when |
|---|---|---|
| [bixverse.plots](https://github.com/GregorLueg/bixverse.plots) | plotting, heavily single cell | any `*_plot_sc()` function: `embedding_plot_sc`, `feature_plot_sc`, `dot_plot_sc`, `stacked_violin_plot_sc`, `violin_plot_sc`, `joint_plot_sc` |
| [bixverse.gpu](https://github.com/GregorLueg/bixverse.gpu) | GPU methods via cubecl/Burn on wgpu, so any GPU | large data. GPU Harmony, PCA and kNN searches. Worth it above a few hundred thousand cells |
| [manifoldsR](https://github.com/GregorLueg/manifoldsR) | UMAP, tSNE, diffusion maps, EVoC clustering | it's already an `Imports`, and all 2D embeddings route through it. `umap_sc()` takes `manifoldsR::params_umap()` |
| [genewalkR](https://github.com/GregorLueg/genewalkR) | gene-gene interaction and regulatory network database | you need a network to diffuse over or to seed a GRN analysis |

`manifoldsR` and `bixverse.plots` are in `Remotes`, so `devtools::install_github`
pulls pinned versions of them automatically. `bixverse.gpu` and `genewalkR` are
separate installs.

## Checking an install

Run the smoke test in `SKILL.md`. It uses synthetic data and needs no network.
If it completes and the printed object shows a populated cache, the Rust side
compiled and linked correctly.

## Keeping this skill current

The skill is versioned with the package. After upgrading bixverse:

```r
bixverse::install_agent_skill(overwrite = TRUE)
```

Other agents work too. `agent = "codex"` (or `"generic"` with an explicit
`dest`) writes the entry file as `AGENTS.md` with the Claude-specific
frontmatter stripped. The `references/` files are plain markdown and portable as
they are. Only Claude Code auto-discovers the skill, so elsewhere you have to
point the agent at the directory yourself.
