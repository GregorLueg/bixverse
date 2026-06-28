# Read a sample of values from a matrix slot in an h5ad file

Handles both sparse (group with a `data` dataset) and dense (direct
dataset) storage. Returns NULL if the slot does not exist.

## Usage

``` r
.read_slot_value_sample(f_path, slot_path, h5_content, n_sample = 10000L)
```
