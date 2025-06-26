library(data.table)
library(igraph)

# Example ontology from your question
test_onto <- data.table(
  parent = c("a", "b", "b", "b", "c"),
  child = c("b", "c", "d", "e", "f")
)

# Function to get all ancestors for each term
get_all_ancestors <- function(ontology_dt) {
  # Create a graph from the ontology
  g <- graph_from_data_frame(ontology_dt, directed = TRUE)

  # Get all unique terms
  all_terms <- unique(c(ontology_dt$parent, ontology_dt$child))

  # For each term, find all its ancestors
  ancestors_list <- list()

  for (term in all_terms) {
    if (term %in% V(g)$name) {
      # Find all nodes that can reach this term (ancestors)
      ancestors <- names(subcomponent(g, term, mode = "in"))
      # Remove the term itself from ancestors
      ancestors <- setdiff(ancestors, term)
      ancestors_list[[term]] <- ancestors
    } else {
      ancestors_list[[term]] <- character(0)
    }
  }

  return(ancestors_list)
}

# Function to calculate S-values for a term's DAG
calculate_s_values <- function(term, ancestors_list, ontology_dt, w = 0.8) {
  # Get ancestors of the term
  ancestors <- ancestors_list[[term]]

  # Create DAG containing term + ancestors
  dag_terms <- c(term, ancestors)

  # Initialize S-values
  s_values <- setNames(rep(0, length(dag_terms)), dag_terms)
  s_values[term] <- 1.0 # The term itself has S-value = 1

  # Build adjacency info for this DAG
  dag_edges <- ontology_dt[parent %in% dag_terms & child %in% dag_terms]

  # Calculate S-values using topological order (from leaves to root)
  # We need to process in reverse topological order
  if (nrow(dag_edges) > 0) {
    g_dag <- graph_from_data_frame(dag_edges, directed = TRUE)
    topo_order <- topo_sort(g_dag, mode = "in") # Reverse topological

    for (node_name in names(topo_order)) {
      if (node_name != term) {
        # Skip the term itself
        # Find children of this node in the DAG
        children <- dag_edges[parent == node_name]$child
        if (length(children) > 0) {
          # S-value is max of w * S-value of children
          s_values[node_name] <- max(w * s_values[children])
        }
      }
    }
  }

  return(s_values)
}

# Function to calculate Wang similarity between two terms
wang_similarity <- function(
  term1,
  term2,
  ancestors_list,
  ontology_dt,
  w = 0.8
) {
  # Calculate S-values for both terms
  s_values_1 <- calculate_s_values(term1, ancestors_list, ontology_dt, w)
  s_values_2 <- calculate_s_values(term2, ancestors_list, ontology_dt, w)

  # Get common terms (intersection of DAGs)
  common_terms <- intersect(names(s_values_1), names(s_values_2))

  # Calculate semantic values (SV) - sum of all S-values
  sv_1 <- sum(s_values_1)
  sv_2 <- sum(s_values_2)

  # Calculate similarity
  if (length(common_terms) == 0) {
    return(0)
  }

  numerator <- sum(s_values_1[common_terms] + s_values_2[common_terms])
  denominator <- sv_1 + sv_2

  similarity <- numerator / denominator

  # Return detailed result
  return(list(
    similarity = similarity,
    common_terms = common_terms,
    s_values_1 = s_values_1,
    s_values_2 = s_values_2,
    sv_1 = sv_1,
    sv_2 = sv_2
  ))
}

# Function to create similarity matrix for all terms
create_similarity_matrix <- function(ontology_dt, w = 0.8) {
  # Get all ancestors
  ancestors_list <- get_all_ancestors(ontology_dt)

  # Get all unique terms
  all_terms <- unique(c(ontology_dt$parent, ontology_dt$child))
  n_terms <- length(all_terms)

  # Initialize similarity matrix
  sim_matrix <- matrix(0, nrow = n_terms, ncol = n_terms)
  rownames(sim_matrix) <- all_terms
  colnames(sim_matrix) <- all_terms

  # Fill the matrix
  for (i in 1:n_terms) {
    for (j in i:n_terms) {
      # Only upper triangle + diagonal
      term1 <- all_terms[i]
      term2 <- all_terms[j]

      if (i == j) {
        sim_matrix[i, j] <- 1.0 # Self-similarity
      } else {
        result <- wang_similarity(term1, term2, ancestors_list, ontology_dt, w)
        sim_matrix[i, j] <- result$similarity
        sim_matrix[j, i] <- result$similarity # Symmetric
      }
    }
  }

  return(sim_matrix)
}

ancestors_list <- get_all_ancestors(test_onto)

result_cd <- wang_similarity("c", "d", ancestors_list, test_onto, w = 0.8)


sim_matrix <- create_similarity_matrix(test_onto, w = 0.8)


test_onto <- data.table::data.table(
  parent = c("a", "b", "b", "b", "c"),
  child = c("b", "c", "d", "e", "f")
)

tictoc::tic()
rs_onto_sim_wang(test_onto$parent, test_onto$child, w = 0.8)
tictoc::toc()
