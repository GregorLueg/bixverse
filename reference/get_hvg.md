# Get the HVG

Returns the HVG indices. Pending class type this are 1-based (for R) or
0-based for Rust.

## Usage

``` r
get_hvg(x)

## S7 method for class <bixverse::MetaCells>
get_hvg(x)

# S3 method for class 'ScMap'
get_hvg(x)

## S7 method for class <bixverse::SingleCells>
get_hvg(x)
```

## Arguments

- x:

  An object to get HVG from.

## Value

Indices of the stored HVG genes.
