use extendr_api::prelude::*;
use std::collections::{BTreeMap, HashMap, HashSet};

use crate::helpers_fgsea::{
    calc_gsea_stats, calc_gsea_stats_helper, calculate_nes_es_pval, GseaResults,
};
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

/// Type alias for intermediary results
/// The first value is the ES score, the second the size and the third
/// the indices of the potential leading edge genes
type GoIntermediaryRes = (f64, usize, Vec<i32>);

/// Return structure of the `process_ontology_level()` ontology function.
#[derive(Clone, Debug)]
pub struct GoElimLevelResults {
    pub go_ids: Vec<String>,
    pub pvals: Vec<f64>,
    pub odds_ratios: Vec<f64>,
    pub hits: Vec<u64>,
    pub gene_set_lengths: Vec<u64>,
}

/// Return structure of the `process_ontology_level()` ontology function.
#[derive(Clone, Debug)]
pub struct GoElimLevelResultsGsea {
    pub go_ids: Vec<String>,
    pub es: Vec<f64>,
    pub nes: Vec<Option<f64>>,
    pub size: Vec<usize>,
    pub pvals: Vec<f64>,
    pub leading_edge: Vec<Vec<i32>>,
}

#[derive(Clone, Debug)]
pub struct GeneOntology<'a> {
    pub go_to_gene: GeneMap,
    pub ancestors: &'a AncestorMap,
    pub levels: &'a LevelMap,
}

impl<'a> GeneOntology<'a> {
    pub fn new(gene_map: GeneMap, ancestor_map: &'a AncestorMap, levels_map: &'a LevelMap) -> Self {
        GeneOntology {
            go_to_gene: gene_map,
            ancestors: ancestor_map,
            levels: levels_map,
        }
    }

    /// Returns the ancestors of a given gene ontology term identifier
    pub fn get_ancestors(&self, id: &String) -> Option<&Vec<String>> {
        self.ancestors.get(id)
    }

    pub fn get_level_ids(&self, id: &String) -> Option<&Vec<String>> {
        self.levels.get(id)
    }

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

/// Structure to hold random permutations for the continuous (fgsea) based
/// enrichment with elimination.
#[derive(Clone, Debug)]
pub struct GeneOntologyRandomPerm<'a> {
    pub random_perm: &'a Vec<Vec<f64>>,
}

impl<'a> GeneOntologyRandomPerm<'a> {
    pub fn new(perm_es: &'a Vec<Vec<f64>>) -> Self {
        GeneOntologyRandomPerm {
            random_perm: perm_es,
        }
    }

    pub fn get_gsea_res_simple<'b>(
        &self,
        pathway_scores: &'b [f64],
        pathway_sizes: &'b [usize],
    ) -> Result<GseaResults<'b>> {
        // Dual lifetimes fun...
        let gsea_batch_res =
            calc_gsea_stats_helper(pathway_scores, pathway_sizes, self.random_perm)?;

        let gsea_res = calculate_nes_es_pval(pathway_scores, pathway_sizes, &gsea_batch_res);

        Ok(gsea_res)
    }
}

///////////////
// Functions //
///////////////

/////////////
// Helpers //
/////////////

/// Take the S7 go_data_class and return the necessary Rust types for further
/// processing.
pub fn prepare_go_data(go_obj: Robj) -> extendr_api::Result<(GeneMap, AncestorMap, LevelMap)> {
    // TODO: Need to do better error handling here... Future me problem
    let go_to_genes = go_obj.get_attrib("go_to_genes").unwrap().as_list().unwrap();
    let ancestors = go_obj.get_attrib("ancestry").unwrap().as_list().unwrap();
    let levels = go_obj.get_attrib("levels").unwrap().as_list().unwrap();

    let go_to_genes = r_list_to_hashmap_set(go_to_genes)?;
    let ancestors = r_list_to_hashmap(ancestors)?;
    let levels = r_list_to_hashmap(levels)?;

    Ok((go_to_genes, ancestors, levels))
}

/////////////////////
// Hypergeom tests //
/////////////////////

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
    let default = vec!["string".to_string()];
    let go_ids = go_obj.get_level_ids(level).unwrap_or(&default);

    let level_data = go_obj.get_genes_list(go_ids);

    let mut level_data_final = HashMap::new();

    for (key, value) in &level_data {
        if value.len() >= min_genes {
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
        if debug {
            println!("This genes are being tested: {:?}", value)
        };
        let gene_set_length = value.len() as u64;
        let hits = target_set.intersection(value).count() as u64;
        if debug {
            println!("Number of hits: {:?}", hits)
        };
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
            "At level {} a total of {} gene ontology terms will be affected by elimination: {:?}",
            level, no_terms, go_to_remove
        );
    }

    for term in go_to_remove.iter() {
        let ancestors = go_obj.get_ancestors(term);
        let ancestors_final: Vec<String> = ancestors.cloned().unwrap_or_else(Vec::new);

        if debug {
            println!(
                "The following ancestors are affected: {:?}",
                ancestors_final
            )
        }

        if let Some(genes_to_remove) = go_obj.get_genes(term) {
            let genes_to_remove = genes_to_remove.clone();
            if debug {
                println!("The following genes will be removed: {:?}", genes_to_remove)
            }
            go_obj.remove_genes(&ancestors_final, &genes_to_remove);
        }

        if debug {
            for ancestor in &ancestors_final {
                let mut binding = HashSet::with_capacity(1);
                binding.insert("no genes left".to_string());
                let new_genes = go_obj.get_genes(ancestor).unwrap_or(&binding);
                if debug {
                    println!("The following genes remain: {:?}", new_genes)
                }
            }
        }
    }

    res
}

/////////////////////
// Continuous test //
/////////////////////

#[allow(clippy::too_many_arguments)]
pub fn process_ontology_level_fgsea_simple(
    stats: &[f64],
    stat_names: &[String],
    gsea_param: f64,
    level: &String,
    go_obj: &mut GeneOntology,
    go_random_perms: &GeneOntologyRandomPerm,
    min_size: usize,
    max_size: usize,
    elim_threshold: f64,
    debug: bool,
) -> Result<GoElimLevelResultsGsea> {
    // Get the identfiers of that level and clean everything up
    let default = vec!["string".to_string()];
    let go_ids = go_obj.get_level_ids(level).unwrap_or(&default);

    // BTreeMap to make sure the order is determistic
    let mut level_data_es: BTreeMap<String, GoIntermediaryRes> = BTreeMap::new();

    for go_id in go_ids {
        if let Some(genes) = go_obj.get_genes(go_id) {
            if genes.len() >= min_size && genes.len() <= max_size {
                // Convert gene names to indices in one step
                let indices: Vec<i32> = stat_names
                    .iter()
                    .enumerate()
                    .filter_map(|(i, s)| {
                        if genes.contains(s) {
                            Some(i as i32)
                        } else {
                            None
                        }
                    })
                    .collect();

                if !indices.is_empty() {
                    let es_res = calc_gsea_stats(stats, &indices, gsea_param, true, false);
                    let size = indices.len();
                    level_data_es.insert(go_id.clone(), (es_res.0, size, es_res.1));
                }
            }
        }
    }

    if debug {
        println!("Generated successful the BTreeMap for level: {:?}", level);
    }

    let mut pathway_scores: Vec<f64> = Vec::with_capacity(level_data_es.len());
    let mut pathway_sizes: Vec<usize> = Vec::with_capacity(level_data_es.len());
    let mut leading_edge_indices = Vec::with_capacity(level_data_es.len());

    for v in level_data_es.values() {
        pathway_scores.push(v.0);
        pathway_sizes.push(v.1);
        leading_edge_indices.push(v.2.clone());
    }

    let level_res: GseaResults<'_> =
        go_random_perms.get_gsea_res_simple(&pathway_scores, &pathway_sizes)?;

    if debug {
        println!(
            "I calculated successfully the random permutations for level: {:?}",
            level
        );
    }

    let go_to_remove: Vec<&String> = level_res
        .pvals
        .iter()
        .zip(level_data_es.keys())
        .filter(|(pval, _)| pval <= &&elim_threshold)
        .map(|(_, go_id)| go_id)
        .collect();

    for term in go_to_remove.iter() {
        let ancestors = go_obj.get_ancestors(term);
        let ancestors_final: Vec<String> = ancestors.cloned().unwrap_or_else(Vec::new);

        if let Some(genes_to_remove) = go_obj.get_genes(term) {
            let genes_to_remove = genes_to_remove.clone();
            go_obj.remove_genes(&ancestors_final, &genes_to_remove);
        }

        if debug {
            for ancestor in &ancestors_final {
                let mut binding = HashSet::with_capacity(1);
                binding.insert("no genes left".to_string());
                let new_genes = go_obj.get_genes(ancestor).unwrap_or(&binding);
                if debug {
                    println!("The following genes remain: {:?}", new_genes)
                }
            }
        }
    }

    Ok(GoElimLevelResultsGsea {
        go_ids: level_data_es.keys().cloned().collect(),
        es: pathway_scores.clone(),
        nes: level_res.nes,
        size: pathway_sizes.clone(),
        pvals: level_res.pvals,
        leading_edge: leading_edge_indices,
    })
}
