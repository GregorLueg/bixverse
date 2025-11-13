# Network diffusion class

This class helps to do diffusion of seed nodes in a single or tied
version over a network, measure the ability of these diffusion vectors
to recall against a gold standard set of nodes and do community
detection within the subset of the network that received the the most
heat from the initial seed genes.

## Usage

``` r
network_diffusions(edge_data_frame, weighted, directed)
```

## Arguments

- edge_data_frame:

  data.table that contains the edge information. It is expected to have
  the columns 'from' and 'to'.

- weighted:

  Boolean. Is the graph weighted. If set to TRUE, the `edge_data_frame`
  needs to have a weight column.

- directed:

  Boolean. Shall the graph be stored as directed.

## Value

Returns the `network_diffusions` class for further operations.

## Properties

- graph:

  igraph. The underlying graph.

- diffusion_res:

  Numeric vector. Contains contains the single or tied diffusion
  results.

- params:

  A (nested) list that will store all the parameters of the applied
  function.

- final_results:

  data.table. Contains final results.
