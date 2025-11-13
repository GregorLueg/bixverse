# Get the embedding from the cache

General wrapper function that can be used to pull out any embedding
stored in the `sc_cache`.

## Usage

``` r
get_embedding(x, embd_name)

# S3 method for class 'sc_cache'
get_embedding(x, embd_name)
```

## Arguments

- x:

  An object to get embedding from

- embd_name:

  String. The name of the embedding to return. The function will throw
  an error if the embedding does not exist.
