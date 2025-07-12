similarity_type = "resnik"
ancestor_list = ancestors
ic_list = information_content

onto_similarities <- rs_onto_semantic_sim_mat(
  sim_type = similarity_type,
  ancestor_list = ancestor_list,
  ic_list = ic_list,
  flat_matrix = FALSE
)

final_sim_mat <- onto_similarities$sim_mat
rownames(final_sim_mat) <- colnames(final_sim_mat) <- onto_similarities$names

final_sim_mat[1:5, 1:5]
