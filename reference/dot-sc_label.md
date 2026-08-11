# Label an artefact for reporting and provenance lookups

Label an artefact for reporting and provenance lookups

## Usage

``` r
.sc_label(modality, artefact, name = NULL)
```

## Arguments

- modality:

  String. The modality the artefact lives in.

- artefact:

  String. One of `c("pca", "embedding", "knn", "snn", "magic")`.

- name:

  String. Embedding name, only used when `artefact = "embedding"`.

## Value

String of the form `"modality:artefact"`, where embeddings use their
name (`"rna:umap"`).
