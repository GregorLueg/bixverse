library(data.table)
library(magrittr)

load("../../workbench/test_communities.rda")

rextendr::document()

rs_rbh_sets(
  module_list = test,
  overlap_coefficient = TRUE,
  min_similarity = 0.2,
  debug = TRUE
)


a <- test$a[c("cluster_23", "cluster_22")]
b <- test$b[c("cluster_1", "cluster_2")]
length(intersect(a, b))
test_rbh = list(a = a, b = b)


map(
  names(a),
  ~ {
    aa = .x
    map(
      names(b),
      ~ {
        bb = .x
        print(paste(
          "names",
          aa,
          bb,
          length(intersect(a[[aa]], b[[bb]]))
        ))
      }
    )
  }
)

a <- test$a[c("cluster_5", "cluster_6", "cluster_7", "cluster_8", "cluster_9")]
b <- test$b[c("cluster_1", "cluster_2")]

rs_rbh_sets(
  module_list = test_rbh,
  overlap_coefficient = TRUE,
  min_similarity = 0.2,
  debug = TRUE
)
