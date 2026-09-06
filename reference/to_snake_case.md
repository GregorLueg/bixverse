# Helper function to transform strings to snake_case

Helper function to transform strings to snake_case

## Usage

``` r
to_snake_case(x, ignore_na = FALSE)
```

## Arguments

- x:

  String. The string (vector) you wish to transform to snake_case.

- ignore_na:

  Boolean. Shall the function just ignore `NA` values and return `NA` at
  this position. Defaults to `FALSE`.

## Value

Returns the string in snake_case format.

## Examples

``` r
# normalise mixed naming conventions
to_snake_case(c("Gene Name", "someRNAValue", "Foo-Bar"))
#> [1] "gene_name"     "some_rnavalue" "foo_bar"      
```
