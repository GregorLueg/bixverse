# Get the parameters that were used.

Extracts parameters from the `BixverseBaseClass` class (or child
classes) and has options to return (pretty) JSONs. This generic also
gets inherited by other classes and can be used to extract parameters.
Also, can dispatch to specific methods for certain S3 classes.

## Usage

``` r
get_params(object, to_json = FALSE, pretty_json = FALSE)

get_params.LdaResult(object, ...)

get_params.LdaKSweepResult(object, ...)

get_params.Hotspot(object, to_json = FALSE, pretty_json = FALSE)

get_params.miloR(object, to_json = FALSE, pretty_json = FALSE)

get_params.ScenicGrn(object, to_json = FALSE, pretty_json = FALSE)

get_params.NmfResult(object, to_json = FALSE, pretty_json = FALSE)

get_params.StabilisedNmfResult(object, to_json = FALSE, pretty_json = FALSE)

get_params.ConsensusNmfResult(object, to_json = FALSE, pretty_json = FALSE)

get_params.DialogueResult(object, ...)

get_params.ScNebula(object, to_json = FALSE, pretty_json = FALSE)
```

## Arguments

- object:

  A class within bixverse that inherits from
  [`BixverseBaseClass()`](https://gregorlueg.github.io/bixverse/reference/BixverseBaseClass.md)
  or defined S3 classes.

- to_json:

  Shall the params be returned as a JSON string.

- pretty_json:

  Shall the params be returned as a pretty JSON string.

- ...:

  Unused, present so the S3 methods sharing this page match the generic.

## Value

Depending on parameters either the R list or a (pretty) JSON string.
