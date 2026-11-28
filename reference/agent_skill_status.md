# Installed bixverse agent skills and their versions

Looks in the user and project locations of every supported agent.
Symlinked installs track their source, so they are left out.

## Usage

``` r
agent_skill_status()
```

## Value

A data.table with `path`, `version` (the stamped package version, `NA`
for an install without a stamp) and `current` (whether it matches the
installed package).

## Examples

``` r
agent_skill_status()
#> Empty data.table (0 rows and 3 cols): path,version,current
```
