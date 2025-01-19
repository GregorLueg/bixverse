use extendr_api::prelude::*;
use std::collections::{HashMap, HashSet};
use crate::utils_r_rust::{r_list_to_hashmap, r_list_to_hashmap_set};
use crate::hypergeom_helpers::*;

// Class, methods
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

  //
  pub fn get_level_ids(
    &self,
    id: &String
  ) -> Option<&Vec<String>> {
    self.levels.get(id)
  }

  // Remove genes from defined sets of genes
  pub fn remove_genes(
    &mut self, 
    ids: &Vec<String>, 
    genes_to_remove: &HashSet<String>
  ) {
    for id in ids.iter() {
      if let Some(gene_set) = self.go_to_gene.get_mut(id) {
          gene_set.retain(|gene| !genes_to_remove.contains(gene));
      }
    }
  }

  // Get the genes for a set of IDs This function will ensure that the keys and ids overlap
  // and return the ids that could be found.
  pub fn get_genes_list(
    &self,
    ids: Vec<String>,
  ) -> (Vec<String>, Vec<&HashSet<String>>) {
    let keys: HashSet<String> = self.go_to_gene
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

  // Get the genes of a specific ID
  pub fn get_genes(
    &self,
    id: &String
  ) -> Option<&HashSet<String>> {
    self.go_to_gene.get(id)
  }
}


// Function to go from R List object to the Gene Ontology structure
pub fn generate_go_structure(
  go_to_genes: List,
  ancestors: List,
  levels: List,
) -> GeneOntology {
  let go_to_genes = r_list_to_hashmap_set(go_to_genes);
  let ancestors = r_list_to_hashmap(ancestors);
  let levels = r_list_to_hashmap(levels);

  let go_obj = GeneOntology {
    go_to_gene: go_to_genes,
    ancestors: ancestors,
    levels: levels,
  };

  return go_obj
}

// Process a given ontology level
pub fn process_ontology_level(
  target_genes: Vec<String>,
  level: &String,
  go_obj: &mut GeneOntology,
  min_genes: i64,
  gene_universe_length: u64,
  elim_threshold: f64,
  debug: bool,
) -> (Vec<String>, Vec<f64>, Vec<f64>, Vec<u64>, Vec<u64>) {
  // Get the identfiers of that level and clean everything up
  let go_ids = go_obj.get_level_ids(&level);
  let go_ids_final: Vec<String> = go_ids
    .map(|v| v.clone())
    .unwrap_or_else(Vec::new);
  let level_data = go_obj.get_genes_list(go_ids_final);
  let mut go_identifiers = level_data.0;
  let mut go_gene_sets = level_data.1;

  // Remove gene ontology terms that do not have a minimum of genes
  let filtered_sets: Vec<(String, &HashSet<String>)> = go_identifiers
    .iter()
    .zip(go_gene_sets.iter())
    .filter(|(_, set)| {
      set.len() as i64 >= min_genes
    })
    .map(|(string, set)| (string.clone(), *set))
    .collect();

  go_identifiers.clear();
  go_gene_sets.clear();

  for (string, set) in filtered_sets {
    go_identifiers.push(string);
    go_gene_sets.push(set);
  }

  // Run the hypergeometric tests on the reduced data

  // Prepare the variables
  let trials = target_genes
    .clone()
    .into_iter()
    .collect::<Vec<_>>()
    .len() as u64;
  let gene_set_lengths = go_gene_sets
    .clone()
    .into_iter()
    .map(|s| {
      s.len() as u64
    })
    .collect::<Vec<u64>>();
  let hits = count_hits_2(go_gene_sets, &target_genes);

  // Calculate p-values and odds ratios
  let pvals: Vec<f64> = hits
    .iter()
    .zip(gene_set_lengths.iter())
    .map(|(hit, gene_set_length)| {
      hypergeom_pval(
        *hit, 
        *gene_set_length, 
        gene_universe_length - *gene_set_length, 
        trials
      )
    })
    .collect();
  let odds_ratios: Vec<f64> = hits
    .iter()
    .zip(gene_set_lengths.iter())
    .map(|(hit, gene_set_length)| {
      hypergeom_odds_ratio(
        *hit,
        *gene_set_length - *hit,
        trials - *hit,
        gene_universe_length - *gene_set_length - trials + *hit,
      )
    })
    .collect();

  // Identify the GO terms were to apply the elimination on (if any)
  let go_to_remove: Vec<String> = go_identifiers
    .iter()
    .zip(pvals.iter())
    .filter(|(_, pval)| {
      pval <= &&elim_threshold
    })
    .map(|(string, _)| string.clone())
    .collect();

  if debug {
    let no_terms = go_to_remove.len();
    println!("At level {} a total of {} gene ontology terms will be affected by elimination.", level, no_terms);
  }

  for term in go_to_remove.iter() {
    let ancestors = go_obj.get_ancestors(term);
    let ancestors_final: Vec<String> = ancestors
      .map(|v| v.clone())
      .unwrap_or_else(Vec::new);
    
    if let Some(genes_to_remove) = go_obj.get_genes(term) {
        let genes_to_remove = genes_to_remove.clone();
        go_obj.remove_genes(&ancestors_final, &genes_to_remove);
    }
  }

  (go_identifiers, pvals, odds_ratios, hits, gene_set_lengths)
}
