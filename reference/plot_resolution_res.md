# Plot the resolution results.

Plots the resolution results (if they can be found in the class). The
x-axis reflects the different resolutions and the y axis the modularity
observed with that resolution.

## Usage

``` r
plot_resolution_res(object, print_head = TRUE, ...)
```

## Arguments

- object:

  The class, either `RbhGraph` or `BulkCoExp`.

- print_head:

  Boolean. Print the Top5 resolution parameters and their meta data.
  Only applicable for `BulkCoExp` objects.

- ...:

  Additional arguments to parse to the functions.

## Value

Plots the result, if the results were found in the class. Otherwise,
throws a warning and returns NULL.
