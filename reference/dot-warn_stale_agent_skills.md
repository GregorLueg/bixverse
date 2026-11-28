# Warn about installed skills that lag the package

Called on attach. A copied skill does not follow a package upgrade, so
an agent would keep reading the API of the old version. Never errors.

## Usage

``` r
.warn_stale_agent_skills()
```

## Value

`NULL`, invisibly.
