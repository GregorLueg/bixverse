# Convert a graph to igraph

Rebuilds the neighbour graph as an undirected igraph with the weights in
the edge attribute `weight` and the build method in the graph attribute
`method`. Costly on large graphs; the package's own consumers work on
the Rust graph directly.

## Usage

``` r
as_igraph(x, ...)

# S3 method for class 'SingleCellConnectivity'
as_igraph(x, ...)

# S3 method for class 'igraph'
as_igraph(x, ...)
```

## Arguments

- x:

  A `SingleCellConnectivity` or an `igraph` (returned as is).

- ...:

  Other parameters.

## Value

An undirected, weighted `igraph`.

## Examples

``` r
conn <- as_sc_connectivity(igraph::make_graph("Zachary"))
as_igraph(conn)
#> IGRAPH 90fcec5 U-W- 34 78 -- 
#> + attr: weight (e/n)
#> + edges from 90fcec5:
#>  [1]  1-- 2  1-- 3  1-- 4  1-- 5  1-- 6  1-- 7  1-- 8  1-- 9  1--11  1--12
#> [11]  1--13  1--14  1--18  1--20  1--22  1--32  2-- 3  2-- 4  2-- 8  2--14
#> [21]  2--18  2--20  2--22  2--31  3-- 4  3-- 8  3-- 9  3--10  3--14  3--28
#> [31]  3--29  3--33  4-- 8  4--13  4--14  5-- 7  5--11  6-- 7  6--11  6--17
#> [41]  7--17  9--31  9--33  9--34 10--34 14--34 15--33 15--34 16--33 16--34
#> [51] 19--33 19--34 20--34 21--33 21--34 23--33 23--34 24--26 24--28 24--30
#> [61] 24--33 24--34 25--26 25--28 25--32 26--32 27--30 27--34 28--34 29--32
#> [71] 29--34 30--33 30--34 31--33 31--34 32--33 32--34 33--34
```
