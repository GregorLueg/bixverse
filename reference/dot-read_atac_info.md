# Read the ATAC state a store records

Rebuilds the kind, count mode and TF-IDF flag from the store's header
and metadata. Stores written before the count mode was recorded are
taken as insertion counts, with a warning.

## Usage

``` r
.read_atac_info(store_dir, lsi = NULL)
```

## Arguments

- store_dir:

  String. The directory holding the store's two files.

- lsi:

  Optional list with the LSI fields of an earlier `AtacInfo`, as saved
  in the memory file.

## Value

An `AtacInfo` object.
