# Libraries ----

library(arrow)
library(data.table)
library(tiledb)
library(tiledbsoma)

# Tests ----

## Can we transform arrow tables directly into data.table? ----

dat <- arrow_table(x = 1:3, y = c("a", "b", "c"))

x <- as.data.table(dat)

## SOMA interactions ----

experiment <- SOMAExperimentOpen("~/Documents/Comp/SOMA/")

experiment$obs$read(coords = 0:9)$concat() |> as.data.table()
