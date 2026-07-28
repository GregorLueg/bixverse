# Guard on detection_method for a BulkCoExp method

Reads `detection_method` from the class params and checks it against the
set of methods the caller is designed to handle. Returns the resolved
method on success, `NULL` on mismatch (with a warning). If
`allow_unset = TRUE`, a `NULL` `detection_method` is treated as a first
invocation and returns `NA_character_` silently.

Caller pattern:


    detection_method <- .assert_bulk_detection_method(object, "correlation-based", "correlation")
    if (is.null(detection_method)) return(object)

## Usage

``` r
.assert_bulk_detection_method(
  object,
  allowed,
  method_label,
  allow_unset = FALSE
)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

- allowed:

  Character vector. Detection-method strings the caller accepts.

- method_label:

  String. Human-readable label used in the warning message (e.g.
  `"correlation"`, `"ICA"`, `"DGRDL"`, `"NMF"`).

- allow_unset:

  Boolean. If `TRUE`, a `NULL` `detection_method` is not an error;
  returns `NA_character_`. Used by finalisers that set
  `detection_method` themselves on first invocation.

## Value

One of:

- the resolved `detection_method` string (proceed)

- `NA_character_` (proceed, first invocation, only when
  `allow_unset = TRUE`)

- `NULL` (mismatch, caller should return object unchanged)
