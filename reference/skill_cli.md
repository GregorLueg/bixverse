# Manage the agent skill from a shell

The command line of the agent skill, the R twin of the Python
`bixverse skill` command:

    Rscript -e 'bixverse::skill_cli()' skill install [--agent claude|codex]
      [--project] [--dest DIR] [--force]
    Rscript -e 'bixverse::skill_cli()' skill status

Run it with the R library the project uses, so the skill matches the
bixverse version the project sees.

## Usage

``` r
skill_cli(args = commandArgs(trailingOnly = TRUE))
```

## Arguments

- args:

  Character. The command line arguments. Defaults to the trailing
  arguments of the R session.

## Value

The exit code, invisibly: `1` from `status` if any install lags the
package, else `0`. Errors on bad arguments.

## Examples

``` r
skill_cli(c("skill", "status"))
#> No bixverse agent skill installed. Run `skill install`.
```
