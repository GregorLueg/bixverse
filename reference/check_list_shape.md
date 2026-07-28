# Check that a value is a list with the required names

Boilerplate guard used at the top of most parameter checkers in this
file: verifies `x` is a list and that all `required_names` are present
in `names(x)`.

## Usage

``` r
check_list_shape(x, required_names)
```

## Arguments

- x:

  The object to check.

- required_names:

  Character vector of names that must be present in `names(x)`.

## Value

`TRUE` if the check was successful, otherwise a checkmate-style error
string.
