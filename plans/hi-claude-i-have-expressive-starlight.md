# Breaking the bixverse <-> bixverse.plots dependency cycle

## Context

bixverse 0.4.7 adds trajectory inference (Palantir, PAGA). bixverse.plots 0.2.3
adds `paga_plot_sc()` plus a `palette` argument across the single cell plots,
and bumps its bound to `bixverse (>= 0.4.7)`. Neither can currently reach main
without the other being there first.

Two couplings cause this:

1. **bixverse's R-CMD-check installs Suggests**, which pulls bixverse.plots
   from an unpinned `Remotes` entry, i.e. whatever `main` is that day. Once
   `Imports: bixverse (>= 0.4.7)` lands on plots' main, bixverse's own check
   can no longer resolve its dependency graph, and 0.4.7 is not tagged yet.
2. **bixverse's pkgdown renders `trajectory_inference.qmd`**, which calls
   `paga_plot_sc()` and `embedding_plot_sc(palette = )`. Neither exists in
   bixverse.plots 0.2.2, and the old `embedding_plot_sc` has no `...`, so the
   palette argument is a hard `unused argument` error rather than a silent drop.

Coupling 1 is the deadlock, and pinning `Remotes` to a tag removes it: you
choose which plots version enters the check graph instead of inheriting
whatever main carries. Note this is not a new dependency, bixverse's check
already installs bixverse.plots from main today and CI is green, so the pin
changes determinism only.

Coupling 2 is not fixable by pinning. It recurs on any release where a vignette
shows a new plot function, and is handled by ordering: the vignette lands in a
docs-only follow-up after both packages are tagged.

Separately, `trajectory_inference.qmd` cannot render on any runner today: it
reads a hardcoded `~/Desktop/zenodo_data/marrow_sample_scseq_counts.h5ad`. The
file is now on Zenodo, so this gets a proper download helper.

Only `trajectory_inference.qmd` breaks against bixverse.plots 0.2.2. The other
seven vignettes using bixverse.plots pass no new arguments and are safe.

## Plan

Three merges, every CI run green. **Nothing here pushes to origin.** Each step
produces local commits on a branch; you drive the pushes and merges, and each
step must be tagged before the next begins.

### Step 1 - bixverse 0.4.7 (site without the trajectory vignette)

Commits on `release-0.4.7`:

**`DESCRIPTION`**
- `Remotes`: pin `GregorLueg/bixverse.plots` to `GregorLueg/bixverse.plots@v0.2.2`.
  Leave `GregorLueg/manifoldsR` alone.
- `Suggests` unchanged, bixverse.plots stays declared.
- Version stays 0.4.7.

v0.2.2 requires only `bixverse (>= 0.4.0)`, so the graph resolves against the
existing v0.4.6 tag with nothing new needed.

**`R/data_single_cell.R`** - add `download_marrow_cd34()` immediately after
`download_cd34_data()` (currently ends line 776). Mirror that function exactly:
`options(timeout = max(300, old_timeout))` with an `on.exit` restore, download
to `tempdir()`, `R.utils::gunzip(dest_file, remove = TRUE)`, return the
uncompressed path.

- url: `https://zenodo.org/records/21892320/files/marrow_sample_scseq_counts.h5ad.gz?download=1`
- dest: `marrow_sample_scseq_counts.h5ad.gz`
- returns: `file.path(temp_dir, "marrow_sample_scseq_counts.h5ad")`

**`vignettes/trajectory_inference.qmd`** (lines 60-68) - replace the hardcoded
path and the `# TODO` comment with `marrow_path <- download_marrow_cd34()`.

**`_pkgdown.yml`**
- Add `download_marrow_cd34` to the download helper block at line 677.
- Temporarily remove `trajectory_inference` from the `articles:` Single Cells
  list (line 124) and its navbar entry (line 71). It goes back in step 3.

Then `devtools::document()` and `air format .`.

Merge to main. `auto-tag.yml` tags `v0.4.7`. R-CMD-check resolves against plots
v0.2.2; pkgdown renders every vignette except the trajectory one.

### Step 2 - bixverse.plots 0.2.3

`DESCRIPTION` only, on `feat-better-feat-plot`:
- Keep the `bixverse (>= 0.4.7)` bump already in the working tree.
- Change `Remotes: GregorLueg/bixverse` to `GregorLueg/bixverse@v0.4.7`.

Its check now installs bixverse from the tag, so the bound is satisfied
honestly rather than by a warning that happens not to fail. Merge to main,
auto-tags `v0.2.3`.

### Step 3 - bixverse docs-only follow-up

Straight to main, no version bump. `auto-tag.yml` fires on the `DESCRIPTION`
change but skips, since `v0.4.7` already exists.

- `DESCRIPTION`: move the `Remotes` pin to `GregorLueg/bixverse.plots@v0.2.3`.
- `_pkgdown.yml`: restore `trajectory_inference` to the `articles:` list and the
  navbar.

pkgdown re-renders with plots 0.2.3 and the trajectory vignette publishes.

### Steady state

Each release bumps the other package's pin to an already-published tag, so
neither ever depends on an untagged version of the other. When a release adds a
vignette that needs a new plot function, that vignette lands in the step 3
follow-up rather than the release merge.

## Critical files

- `/Users/gregorlueg/repos/shared/bixverse/DESCRIPTION`
- `/Users/gregorlueg/repos/shared/bixverse/R/data_single_cell.R` (~line 776,
  next to `download_cd34_data`)
- `/Users/gregorlueg/repos/shared/bixverse/vignettes/trajectory_inference.qmd`
  (lines 60-68)
- `/Users/gregorlueg/repos/shared/bixverse/_pkgdown.yml` (lines 71, 124, 677)
- `/Users/gregorlueg/repos/shared/bixverse.plots/DESCRIPTION`

## Verification

Locally, before step 1:

```r
devtools::document()
devtools::install()
```

```sh
air format .
```

Confirm the helper fetches and gunzips:

```r
p <- bixverse::download_marrow_cd34()
bixverse::read_h5ad_metadata(p)$dims
```

Render the trajectory vignette against your local bixverse.plots checkout
before trusting CI, since this is the only thing step 3 turns on:

```r
quarto::quarto_render("vignettes/trajectory_inference.qmd")
```

After each merge, on GitHub:

- Step 1: R-CMD-check green on both runners, pkgdown green, `v0.4.7` tagged.
  Check the dependency-install log shows bixverse.plots 0.2.2, not main.
- Step 2: bixverse.plots R-CMD-check green with bixverse resolved from the
  `v0.4.7` tag, `v0.2.3` tagged.
- Step 3: pkgdown green and `articles/trajectory_inference.html` live with PAGA
  plots rendered.

## Notes, not in scope

- bixverse.plots calls `bixverse:::per_cell_qc_outlier` and
  `bixverse:::prep_stats_pathways`. Triple-colon coupling to internals breaks
  silently on a refactor and no version bound catches it.
- The trajectory vignette now downloads a large h5ad on every pkgdown run.
  Worth watching the job duration.
- `auto-tag.yml` in both repos references `steps.release_notes.outputs.notes`
  with no such step defined, so release bodies are always empty, and
  `actions/create-release@v1` is archived.
- bixverse.plots' `R-cmd-check.yml` has the `Setup Rust` step duplicated
  verbatim.
- Neither pkgdown workflow sets up the quarto CLI despite
  `VignetteBuilder: quarto`. Green today, so leave it alone.
