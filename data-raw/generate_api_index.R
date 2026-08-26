# generate the agent skill api index -------------------------------------------
#
# Emits inst/skills/bixverse/references/api-index.md from _pkgdown.yml (the
# taxonomy) and man/*.Rd (the titles). Run from the package root:
#
#   Rscript data-raw/generate_api_index.R
#
# Re-run on every release, otherwise the index drifts from the code.

## libraries -------------------------------------------------------------------

library(yaml)
library(tools)

## config ----------------------------------------------------------------------

out_path <- file.path(
  "inst",
  "skills",
  "bixverse",
  "references",
  "api-index.md"
)

# these get one line instead of ~92, an agent should be calling the R wrapper
collapse_sections <- "Rust wrappers"

## helpers ---------------------------------------------------------------------

#' Pull the character content of the first Rd block with a given tag
rd_tag <- function(rd, tag) {
  hits <- Filter(\(x) identical(attr(x, "Rd_tag"), tag), rd)
  if (length(hits) == 0) {
    return(NA_character_)
  }
  squish(paste(unlist(hits[[1]]), collapse = ""))
}

#' Pull every value for a tag that can repeat (aliases)
rd_tags_all <- function(rd, tag) {
  hits <- Filter(\(x) identical(attr(x, "Rd_tag"), tag), rd)
  vapply(hits, \(x) squish(paste(unlist(x), collapse = "")), character(1))
}

#' Collapse whitespace and strip the markup that survives unlist()
squish <- function(x) {
  x <- gsub("\\\\[a-zA-Z]+\\{|\\}", "", x)
  trimws(gsub("[[:space:]]+", " ", x))
}

## rd titles -------------------------------------------------------------------

rd_files <- list.files("man", pattern = "\\.Rd$", full.names = TRUE)
stopifnot(
  "No Rd files found, run devtools::document() first" = length(
    rd_files
  ) >
    0
)

alias_to_title <- new.env(parent = emptyenv())

for (f in rd_files) {
  rd <- tryCatch(tools::parse_Rd(f), error = \(e) NULL)
  if (is.null(rd)) {
    warning(sprintf("Could not parse %s", f))
    next
  }
  title <- rd_tag(rd, "\\title")
  aliases <- c(rd_tag(rd, "\\name"), rd_tags_all(rd, "\\alias"))
  for (a in unique(aliases[!is.na(aliases)])) {
    assign(a, title, envir = alias_to_title)
  }
}

lookup_title <- function(topic) {
  if (exists(topic, envir = alias_to_title, inherits = FALSE)) {
    get(topic, envir = alias_to_title)
  } else {
    NA_character_
  }
}

## exports ---------------------------------------------------------------------

ns <- readLines("NAMESPACE")
exports <- grep("^export\\(", ns, value = TRUE)
exports <- gsub("^export\\(|\\)$", "", exports)
exports <- gsub('^"|"$|^`|`$', "", exports)

## pkgdown taxonomy ------------------------------------------------------------

pkgdown <- yaml::read_yaml("_pkgdown.yml")
sections <- pkgdown$reference

indexed <- unlist(lapply(sections, \(s) s$contents), use.names = FALSE)

missing_from_pkgdown <- setdiff(exports, indexed)
missing_from_ns <- setdiff(indexed, exports)

if (length(missing_from_ns) > 0) {
  warning(sprintf(
    "In _pkgdown.yml but not exported: %s",
    paste(missing_from_ns, collapse = ", ")
  ))
}

## emit ------------------------------------------------------------------------

lines <- c(
  "# bixverse API index",
  "",
  paste(
    "Every documented entry point, grouped the way the package website groups",
    "them. Generated from `_pkgdown.yml` and `man/*.Rd` by",
    "`data-raw/generate_api_index.R`. Do not edit by hand."
  ),
  "",
  "Grep this file to check whether a function exists before calling it.",
  ""
)

for (s in sections) {
  title <- s$title
  contents <- s$contents
  if (is.null(contents)) {
    next
  }

  lines <- c(lines, sprintf("## %s", title), "")

  if (!is.null(s$desc)) {
    lines <- c(lines, squish(s$desc), "")
  }

  if (title %in% collapse_sections) {
    lines <- c(
      lines,
      sprintf(
        paste(
          "%d `rs_*` functions are exposed here. They are the raw extendr",
          "bindings with no input validation. Use the R wrapper instead; only",
          "reach for these if you are building on top of bixverse and know",
          "exactly what you are doing."
        ),
        length(contents)
      ),
      ""
    )
    next
  }

  for (topic in contents) {
    tl <- lookup_title(topic)
    lines <- c(
      lines,
      if (is.na(tl)) {
        sprintf("- `%s`", topic)
      } else {
        sprintf("- `%s`: %s", topic, tl)
      }
    )
  }
  lines <- c(lines, "")
}

if (length(missing_from_pkgdown) > 0) {
  lines <- c(
    lines,
    "## Not on the package website",
    "",
    paste(
      "Exported but absent from `_pkgdown.yml`. Mostly internal constructors",
      "and `rs_*` bindings that take on-disk streaming input. Usable, but",
      "undocumented on the website, so read the roxygen with `?fn` first."
    ),
    ""
  )
  for (topic in sort(missing_from_pkgdown)) {
    tl <- lookup_title(topic)
    lines <- c(
      lines,
      if (is.na(tl)) {
        sprintf("- `%s`", topic)
      } else {
        sprintf("- `%s`: %s", topic, tl)
      }
    )
  }
  lines <- c(lines, "")
}

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
writeLines(lines, out_path)

message(sprintf(
  "Wrote %s: %d sections, %d indexed topics, %d unindexed exports",
  out_path,
  length(sections),
  length(indexed),
  length(missing_from_pkgdown)
))
