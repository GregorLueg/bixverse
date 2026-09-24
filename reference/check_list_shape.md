# Check that a value is a list with the required names

Boilerplate guard at the top of every generated parameter checker:
verifies `x` is a list and that `required_names` are present in
`names(x)`.

## Usage

``` r
check_list_shape(x, required_names, strict = FALSE)
```

## Arguments

- x:

  The object to check.

- required_names:

  Character vector of names that must be present in `names(x)`.

- strict:

  Boolean. `TRUE` additionally rejects names that are not in
  `required_names`. Defaults to `FALSE`.

## Value

`TRUE` if the check was successful, otherwise a checkmate-style error
string.
