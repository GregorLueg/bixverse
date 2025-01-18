use extendr_api::prelude::*;
use std::collections::{HashMap, HashSet};
use crate::utils_r_rust::{r_list_to_hashmap, r_list_to_hashmap_set};

// Class, methods and generator
pub struct GeneOntology {
    pub go_to_gene: HashMap<String, HashSet<String>>,
    pub ancestors: HashMap<String, Vec<String>>,
    pub levels: HashMap<String, Vec<String>>,
}

impl GeneOntology {
  // Get the ancestors for a given id
  pub fn get_ancestors(
    &self, 
    id: &String
  ) -> Option<&Vec<String>> {
    self.ancestors.get(id)
  }

  // Remove genes from defined sets of genes
  pub fn remove_genes(
    &mut self, 
    ids: Vec<String>, 
    genes_to_remove: HashSet<String>
  ) {
    for id in ids.iter() {
      if let Some(gene_set) = self.go_to_gene.get_mut(id) {
          gene_set.retain(|gene| !genes_to_remove.contains(gene));
      }
    }
  }

  // Get the genes for a set of IDs
  pub fn get_genes(
    &self,
    ids: Vec<String>,
  ) -> (Vec<String>, Vec<&HashSet<String>>) {
    let keys: HashSet<_> = self.go_to_gene
      .clone()
      .into_keys()
      .collect();

    let id_keys: HashSet<_> = ids.into_iter().collect();

    let ids_final: Vec<String> = id_keys
      .intersection(&keys)
      .into_iter()
      .map(|s| s.to_string())
      .collect();

    let gene_sets: Vec<_> = ids_final
      .iter()
      .filter_map(|s| {
        self.go_to_gene.get(s)
      })
      .collect();

    (ids_final, gene_sets)
  }
}

pub fn generate_go_structure(
  go_to_genes: List,
  ancestors: List,
  levels: List,
) -> GeneOntology {
  let go_to_genes = r_list_to_hashmap_set(go_to_genes);
  let ancestors = r_list_to_hashmap(ancestors);
  let levels = r_list_to_hashmap(levels);

  let mut go_obj = GeneOntology {
        go_to_gene: go_to_genes,
        ancestors: ancestors,
        levels: levels,
  };

  return go_obj
}

extendr_module! {
    mod geom_elim_fun;
}