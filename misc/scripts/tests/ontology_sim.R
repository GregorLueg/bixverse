hpo$ancestors[['HP:0000001']]

ontology <- hpo


tab <- table(factor(
  unlist(use.names = FALSE, ontology$ancestors),
  levels = ontology$id
))

ic <- setNames(nm = ontology$id, -log(as.integer(tab) / length(ontology$id)))

represented <- purrr::map_lgl(names(ontology$ancestors), \(x) {
  x %in% names(information_content)
})

names(ontology$ancestors)[!represented]

hpo$parents['HP:0000057']

hpo$ancestors['HP:0000057']

hpo$obsolete['HP:0000057']

which[names(ontology$ancestors) %in% names(information_content)]


information_content['HP:0000001']

ontoSimiliarity_information_content['HP:0000001']

ancestor_list <- ancestors

ontology$ancestors[['HP:0000001']]

ancestors[['HP:0000001']]

ancestor_list <- descendants

tab_2 <- table(factor(
  unlist(use.names = FALSE, ancestors),
  levels = names(ancestors)
))

tab_3 <- purrr::map_dbl(ancestor_list, length)

information_content_v3 <- stats::setNames(
  -log(as.integer(tab) / length(ancestor_list)),
  nm = names(ancestor_list)
)

tab_2['HP:0000001']

tab_3['HP:0000001']

tab['HP:0000001']

ic['HP:0000002']

tab['HP:0000001']
