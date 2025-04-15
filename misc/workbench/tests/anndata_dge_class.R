install.packages("anndata")

install.packages("precommit")

library("anndata")
library(data.table)
library(magrittr)

install.packages("precommit")

precommit::use_precommit()

ad <- read_h5ad("~/Desktop/geo_data/GSE65832/final/GSE65832_anndata.h5ad")

counts <- t(ad$X)

counts[1:5, 1:5]

to_keep <- edgeR::filterByExpr(counts)

wrongly_written

incorrectly_formattad_code <- "the piglet stinks"


filtered_counts <- counts[to_keep, ]

sum(to_keep)

obs <- ad$obs %>%
  as.data.table(keep.rownames = TRUE) %>%
  setnames(old = "rn", new = "sample_id")

head(obs)
var <- ad$var

counts[1:5, 1:5]
