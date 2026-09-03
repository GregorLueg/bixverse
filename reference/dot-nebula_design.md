# Build the NEBULA design from an obs table

Evaluates the design formula against the obs table of the selected cells
(or meta cells) and drops anything with a missing design or subject
value. The subject labels come back as a factor so the caller can map
them onto the full store.

## Usage

``` r
.nebula_design(obs, design, subject_col)
```

## Arguments

- obs:

  data.table. The obs table of the selected cells, holding at least the
  subject column and every variable in `design`.

- design:

  Formula. The experimental design, evaluated against `obs`.

- subject_col:

  String. The column of `obs` holding the subject identifier, i.e. what
  the random effect is over.

## Value

A list with `obs` (rows actually used), `design_mat` and `subject_fct`.
