# Strip a YAML frontmatter block from markdown lines

The Claude skill format puts `name` and `description` in frontmatter.
Other agents do not read it, so it becomes noise at the top of the file.

## Usage

``` r
.strip_frontmatter(lines)
```

## Arguments

- lines:

  Character. The file contents.

## Value

Character. The contents without a leading frontmatter block.
