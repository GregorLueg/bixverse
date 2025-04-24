use extendr_api::prelude::*;
use std::collections::{HashMap, HashSet};

use crate::helpers_hypergeom::*;
use crate::utils_r_rust::{r_list_to_hashmap, r_list_to_hashmap_set};

///////////////////////
// Types & Structure //
///////////////////////

/// Type alias for the go identifier to gene Hashmap
type GeneMap = HashMap<String, HashSet<String>>;

/// Type alias for the ancestor to go identifier HashMap
type AncestorMap = HashMap<String, Vec<String>>;

/// Type alias for the ontology level to go identifier HashMap
type LevelMap = HashMap<String, Vec<String>>;

/// Return structure of the `process_ontology_level()` ontology function.
pub struct GoElimLevelResults {
    pub go_ids: Vec<String>,
    pub pvals: Vec<f64>,
    pub odds_ratios: Vec<f64>,
    pub hits: Vec<u64>,
    pub gene_set_lengths: Vec<u64>,
}

/// Structure that contains the gene ontology and key functions to do apply the
/// elimination method.
pub struct GeneOntology {
    pub go_to_gene: GeneMap,
    pub ancestors: AncestorMap,
    pub levels: LevelMap,
}

impl GeneOntology {
    /// Returns the ancestors of a given gene ontology term identifier
    pub fn get_ancestors(&self, id: &String) -> Option<&Vec<String>> {
        self.ancestors.get(id)
    }

    /// Returns the gene ontology term identifiers for a given level of the
    /// ontology.
    pub fn get_level_ids(&self, id: &String) -> Option<&Vec<String>> {
        self.levels.get(id)
    }

    /// Remove genes from defined sets of genes
    pub fn remove_genes(&mut self, ids: &[String], genes_to_remove: &HashSet<String>) {
        for id in ids.iter() {
            if let Some(gene_set) = self.go_to_gene.get_mut(id) {
                gene_set.retain(|gene| !genes_to_remove.contains(gene));
            }
        }
    }

    /// Get the genes based on an array of Strings.
    pub fn get_genes_list(&self, ids: &[String]) -> HashMap<String, &HashSet<String>> {
        let mut to_ret = HashMap::new();

        for id in ids.iter() {
            if self.go_to_gene.contains_key(id) {
                to_ret.insert(id.to_string(), self.go_to_gene.get(id).unwrap());
            }
        }

        to_ret
    }

    /// Get the genes for one specific ID
    pub fn get_genes(&self, id: &String) -> Option<&HashSet<String>> {
        self.go_to_gene.get(id)
    }
}

///////////////
// Functions //
///////////////

/// Take the S7 go_data_class and return the necessary Rust types for further
/// processing.
pub fn prepare_go_data(go_obj: Robj) -> extendr_api::Result<(GeneMap, AncestorMap, LevelMap)> {
    let go_to_genes = go_obj.get_attrib("go_to_genes").unwrap().as_list().unwrap();
    let ancestors = go_obj.get_attrib("ancestry").unwrap().as_list().unwrap();
    let levels = go_obj.get_attrib("levels").unwrap().as_list().unwrap();

    let go_to_genes = r_list_to_hashmap_set(go_to_genes)?;
    let ancestors = r_list_to_hashmap(ancestors)?;
    let levels = r_list_to_hashmap(levels)?;

    Ok((go_to_genes, ancestors, levels))
}

/// Process a given ontology level
pub fn process_ontology_level(
    target_genes: &[String],
    level: &String,
    go_obj: &mut GeneOntology,
    min_genes: usize,
    gene_universe_length: u64,
    elim_threshold: f64,
    debug: bool,
) -> GoElimLevelResults {
    // Get the identfiers of that level and clean everything up
    let go_ids = go_obj.get_level_ids(level);
    let go_ids_final: &Vec<String> = go_ids.unwrap();
    let level_data = go_obj.get_genes_list(go_ids_final);

    let mut level_data_final = HashMap::new();

    for (key, value) in &level_data {
        if value.len() < min_genes {
            level_data_final.insert(key.clone(), value);
        }
    }

    let trials = target_genes.len() as u64;
    let mut target_set = HashSet::new();
    for s in target_genes {
        target_set.insert(s.clone());
    }

    let size = level_data_final.len();

    let mut go_ids = Vec::with_capacity(size);
    let mut hits_vec = Vec::with_capacity(size);
    let mut pvals = Vec::with_capacity(size);
    let mut odds_ratios = Vec::with_capacity(size);
    let mut gene_set_lengths = Vec::with_capacity(size);

    for (key, value) in level_data_final {
        let gene_set_length = value.len() as u64;
        let hits = target_set.intersection(value).count() as u64;
        let q = hits as i64 - 1;
        let pval = if q > 0 {
            hypergeom_pval(
                q as u64,
                gene_set_length,
                gene_universe_length - gene_set_length,
                trials,
            )
        } else {
            1.0
        };
        let odds_ratio = hypergeom_odds_ratio(
            hits,
            gene_set_length - hits,
            trials - hits,
            gene_universe_length - gene_set_length - trials + hits,
        );
        go_ids.push(key.clone());
        hits_vec.push(hits);
        pvals.push(pval);
        odds_ratios.push(odds_ratio);
        gene_set_lengths.push(gene_set_length);
    }

    let res = GoElimLevelResults {
        go_ids,
        pvals,
        odds_ratios,
        hits: hits_vec,
        gene_set_lengths,
    };

    // Identify the GO terms were to apply the elimination on (if any)
    let go_to_remove: &Vec<_> = &res
        .go_ids
        .iter()
        .zip(&res.pvals)
        .filter(|(_, pval)| pval <= &&elim_threshold)
        .map(|(string, _)| string.clone())
        .collect();

    if debug {
        let no_terms = go_to_remove.len();
        println!(
            "At level {} a total of {} gene ontology terms will be affected by elimination.",
            level, no_terms
        );
    }

    for term in go_to_remove.iter() {
        let ancestors = go_obj.get_ancestors(term);
        let ancestors_final: Vec<String> = ancestors.cloned().unwrap_or_else(Vec::new);

        if let Some(genes_to_remove) = go_obj.get_genes(term) {
            let genes_to_remove = genes_to_remove.clone();
            go_obj.remove_genes(&ancestors_final, &genes_to_remove);
        }
    }

    res
}
