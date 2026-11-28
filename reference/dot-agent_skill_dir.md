# The skills directory an agent discovers skills in

The skills directory an agent discovers skills in

## Usage

``` r
.agent_skill_dir(agent, scope)
```

## Arguments

- agent:

  String. One of `c("claude", "codex")`.

- scope:

  String. One of `c("user", "project")`.

## Value

String. The directory, without the `bixverse` subdirectory. Project
paths are relative to the working directory.
