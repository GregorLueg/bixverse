# Lay out an existing Bonsai tree

**\[experimental\]** Computes a new 2D layout for a finished tree
without searching again. The tree is renumbered internally, so inferred
ancestors can come back with different indices than they went in with;
leaves keep theirs. Replace the whole node table with the output.

## Usage

``` r
rs_bonsai_layout(parent, branch, n_leaves, layout, hyperbolic)
```

## Arguments

- parent:

  Integer. Parent of each node (0-indexed!), negative for the root.

- branch:

  Numeric. Branch length above each node.

- n_leaves:

  Integer. Number of leaves, which occupy the first `n_leaves` nodes.

- layout:

  String. One of `c("equal_angle", "equal_daylight", "dendrogram")`.

- hyperbolic:

  Boolean. Project onto the hyperbolic disk.

## Value

A list with `parent` (0-indexed!, `-1` for the root), `branch`, `x` and
`y`, all indexed by the tree's own node numbering.
