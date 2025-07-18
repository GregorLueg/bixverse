use extendr_api::prelude::*;

use rustc_hash::{FxBuildHasher, FxHashMap, FxHashSet};
use std::collections::BTreeMap;

use crate::helpers::fgsea::{
    calc_gsea_stats, calc_gsea_stats_wrapper, calculate_nes_es_pval, GseaBatchResults, GseaParams,
    GseaResults,
};
use crate::helpers::hypergeom::*;
use crate::utils::general::flatten_vector;
use crate::utils::r_rust_interface::{r_list_to_hashmap, r_list_to_hashmap_set};
use crate::utils::utils_stats::calc_fdr;

///////////////////////
// Types & Structure //
///////////////////////

/// Type alias for the go identifier to gene Hashmap
type GeneMap = FxHashMap<String, FxHashSet<String>>;

/// Type alias for the ancestor to go identifier HashMap
type AncestorMap = FxHashMap<String, Vec<String>>;

/// Type alias for the ontology level to go identifier HashMap
type LevelMap = FxHashMap<String, Vec<String>>;

/// Type alias for intermediary results
///
/// ### Fields
///
/// * `0` - ES scores
/// * `1` - Size
/// * `2` - Indices of the leading edge genes.
type GoIntermediaryRes = (f64, usize, Vec<i32>);

/// Return structure of the `process_ontology_level()` ontology function.
///
/// ### Fields
///
/// * `go_ids` - GO term identifiers.
/// * `pvals` - p-values for the terms.
/// * `odds_ratio` - Calculated odds ratios for the terms.
/// * `hits` - Number of intersecting genes with the target gene set
/// * `gene_set_lengths` - Length of the gene set after elimination.
#[derive(Clone, Debug)]
pub struct GoElimLevelResults {
    pub go_ids: Vec<String>,
    pub pvals: Vec<f64>,
    pub odds_ratios: Vec<f64>,
    pub hits: Vec<usize>,
    pub gene_set_lengths: Vec<usize>,
}

/// Final return structure after filtering
///
/// ### Fields
///
/// * `go_ids` - GO term identifiers.
/// * `pvals` - p-values for the terms.
/// * `fdr` - the FDRs.
/// * `odds_ratio` - Calculated odds ratios for the terms.
/// * `hits` - Number of intersecting genes with the target gene set
/// * `gene_set_lengths` - Length of the gene set after elimination.
#[derive(Clone, Debug)]
pub struct GoElimFinalResults {
    pub go_ids: Vec<String>,
    pub pvals: Vec<f64>,
    pub fdr: Vec<f64>,
    pub odds_ratios: Vec<f64>,
    pub hits: Vec<usize>,
    pub gs_length: Vec<usize>,
}

/// Return structure of the `process_ontology_level_fgsea_simple()` ontology function.
///
/// ### Fields
///
/// * `go_ids` - GO term identifiers.
/// * `es` - Enrichment Scores for the terms.
/// * `nes` - Normalised enrichment Scores for the terms.
/// * `size` - The sizes of the terms.
/// * `pvals` - The p-value for the terms.
/// * `n_more_extreme` - The number of permuted values that were larger (or smaller).
/// * `ge_zero` - The number of permuted values that were ≥ 0.
/// * `le_zero` - The number of permuted values that were ≤ 0.
/// * `leading_edge` - Indices of the leading edge genes of the term.
#[derive(Clone, Debug)]
pub struct GoElimLevelResultsGsea {
    pub go_ids: Vec<String>,
    pub es: Vec<f64>,
    pub nes: Vec<Option<f64>>,
    pub size: Vec<usize>,
    pub pvals: Vec<f64>,
    pub n_more_extreme: Vec<usize>,
    pub ge_zero: Vec<usize>,
    pub le_zero: Vec<usize>,
    pub leading_edge: Vec<Vec<i32>>,
}

///////////////////////
// GO data structure //
///////////////////////

/// Structure for the GeneOntology data
///
/// ### Fields
///
/// * `go_to_gene` - A HashMap with the term to gene associations
/// * `ancestors` - A (borrowed) HashMap with the term to ancestor associations.
/// * `levels` - A (borrowed) HashMap with the ontology level to term associations.
#[derive(Clone, Debug)]
pub struct GeneOntology<'a> {
    pub go_to_gene: GeneMap,
    pub ancestors: &'a AncestorMap,
    pub levels: &'a LevelMap,
}

impl<'a> GeneOntology<'a> {
    /// Create a new GeneOntology Structure
    ///
    /// ### Params
    ///
    /// * `gene_map` - The HashMap with the term to gene associations.
    /// * `ancestors` - The HashMap with the term to ancestor associations.
    /// * `levels` - The HashMap with the level to term associations.
    ///
    /// ### Returns
    ///
    /// Initialised structure
    pub fn new(gene_map: GeneMap, ancestor_map: &'a AncestorMap, levels_map: &'a LevelMap) -> Self {
        GeneOntology {
            go_to_gene: gene_map,
            ancestors: ancestor_map,
            levels: levels_map,
        }
    }

    /// Returns the ancestors of a given gene ontology term identifier
    ///
    /// ### Params
    ///
    /// * `id` - The identifier of the term for which to retrieve the ancestors
    ///
    /// ### Returns
    ///
    /// The ancestors of that term.
    pub fn get_ancestors(&self, id: &String) -> Option<&Vec<String>> {
        self.ancestors.get(id)
    }

    /// Returns the terms for a given level identifier
    ///
    /// ### Params
    ///
    /// * `id` - The level identifier for which to return the GO terms.
    ///
    /// ### Returns
    ///
    /// The GO terms at this level.
    pub fn get_level_ids(&self, id: &String) -> Option<&Vec<String>> {
        self.levels.get(id)
    }

    /// Eliminate genes in a subset of GO terms
    ///
    /// This function implements the elimination logic. For a set of provided
    /// terms it will remove the genes (if they are part of the term) for the
    /// provided GO terms.
    ///
    /// ### Params
    ///
    /// * `ids` - The GO terms in which to remove the genes.
    /// * `genes_to_remove` - HashSet of the genes to remove.
    pub fn remove_genes(&mut self, ids: &[String], genes_to_remove: &FxHashSet<String>) {
        for id in ids.iter() {
            if let Some(gene_set) = self.go_to_gene.get_mut(id) {
                gene_set.retain(|gene| !genes_to_remove.contains(gene));
            }
        }
    }

    /// Get the term to gene associations
    ///
    /// ### Param
    ///
    /// * `ids` - The GO terms for which to return the genes
    ///
    /// ### Returns
    ///
    /// A HashMap with the genes as HashSets.
    pub fn get_genes_list(&self, ids: &[String]) -> FxHashMap<String, &FxHashSet<String>> {
        let mut to_ret = FxHashMap::with_capacity_and_hasher(ids.len(), FxBuildHasher);

        for id in ids.iter() {
            if let Some(gene_set) = self.go_to_gene.get(id) {
                to_ret.insert(id.clone(), gene_set);
            }
        }

        to_ret
    }

    /// Return the genes for a specific term
    ///
    /// ### Params
    ///
    /// * `id` - The GO term for which to return the genes
    ///
    /// ### Returns
    ///
    /// The HashSet with the genes.
    pub fn get_genes(&self, id: &String) -> Option<&FxHashSet<String>> {
        self.go_to_gene.get(id)
    }
}

//////////////////////////////
// GO permutation structure //
//////////////////////////////

/// Structure to hold random permutations for the continuous (fgsea)
///
/// ### Fields
///
/// * `random_perm` - A borrowed vector of vectors with the permutation ES.
#[derive(Clone, Debug)]
pub struct GeneOntologyRandomPerm<'a> {
    pub random_perm: &'a Vec<Vec<f64>>,
}

impl<'a> GeneOntologyRandomPerm<'a> {
    /// Initialise the structure
    ///
    /// ### Params
    ///
    /// * `perm_es` - A vector of vectors with the permutation results
    ///
    /// ### Returns
    ///
    /// Initialised structure
    pub fn new(perm_es: &'a Vec<Vec<f64>>) -> Self {
        GeneOntologyRandomPerm {
            random_perm: perm_es,
        }
    }

    /// Return GseaResults (simple)
    ///
    /// ### Params
    ///
    /// * `pathway_scores` - The pathway scores for which to estimate the permutation
    ///                      based stats
    /// * `pathway_sizes` - The corresponding pathway sizes.
    ///
    /// ### Returns
    ///
    /// A `GseaResults` structure with all of the statistics.
    pub fn get_gsea_res_simple<'b>(
        &self,
        pathway_scores: &'b [f64],
        pathway_sizes: &'b [usize],
    ) -> Result<GseaResults<'b>> {
        // Dual lifetimes fun...
        let gsea_batch_res: GseaBatchResults =
            calc_gsea_stats_wrapper(pathway_scores, pathway_sizes, self.random_perm);

        let gsea_res: GseaResults<'_> =
            calculate_nes_es_pval(pathway_scores, pathway_sizes, &gsea_batch_res);

        Ok(gsea_res)
    }
}

///////////////
// Functions //
///////////////

/////////////
// Helpers //
/////////////

/// Transform the S7 class into the needed data
///
/// ### Params
///
/// * `go_obj` - The S7 class storing the Gene Ontology data for the elimination methods
///
/// ### Returns
///
/// A tuple with the `GeneMap`, `AncestorMap` and `LevelMap`.
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

/// Finalise the GO enrichment results from a hypergeometric test
///
/// ### Params
///
/// * `go_res` - Slice of `GoElimLevelResults` structures
/// * `min_overlap` - Optional minimum overlap.
/// * `fdr_threshold` - Optional fdr threshold.
///
/// ### Returns
///
/// A `GoElimFinalResults` results structure
pub fn finalise_go_res(
    go_res: &[GoElimLevelResults],
    min_overlap: Option<usize>,
    fdr_threshold: Option<f64>,
) -> GoElimFinalResults {
    let n = go_res.iter().map(|x| x.go_ids.len()).sum::<usize>();

    let mut go_ids: Vec<Vec<String>> = Vec::with_capacity(n);
    let mut pvals: Vec<Vec<f64>> = Vec::with_capacity(n);
    let mut hits: Vec<Vec<usize>> = Vec::with_capacity(n);
    let mut odds_ratios: Vec<Vec<f64>> = Vec::with_capacity(n);
    let mut gs_lengths: Vec<Vec<usize>> = Vec::with_capacity(n);

    for res in go_res {
        go_ids.push(res.go_ids.clone());
        pvals.push(res.pvals.clone());
        hits.push(res.hits.clone());
        odds_ratios.push(res.odds_ratios.clone());
        gs_lengths.push(res.gene_set_lengths.clone());
    }

    let go_ids = flatten_vector(go_ids);
    let pvals = flatten_vector(pvals);
    let fdr = calc_fdr(&pvals);
    let hits = flatten_vector(hits);
    let odds_ratios = flatten_vector(odds_ratios);
    let gs_lengths = flatten_vector(gs_lengths);

    let to_keep: Vec<usize> = (0..n)
        .filter(|i| {
            if let Some(min_overlap) = min_overlap {
                if hits[*i] < min_overlap {
                    return false;
                }
            }
            if let Some(fdr_threshold) = fdr_threshold {
                if fdr[*i] > fdr_threshold {
                    return false;
                }
            }
            true
        })
        .collect();

    GoElimFinalResults {
        go_ids: to_keep.iter().map(|i| go_ids[*i].clone()).collect(),
        pvals: to_keep.iter().map(|i| pvals[*i]).collect(),
        fdr: to_keep.iter().map(|i| fdr[*i]).collect(),
        odds_ratios: to_keep.iter().map(|i| odds_ratios[*i]).collect(),
        hits: to_keep.iter().map(|i| hits[*i]).collect(),
        gs_length: to_keep.iter().map(|i| gs_lengths[*i]).collect(),
    }
}

/////////////////////
// Hypergeom tests //
/////////////////////

/// Process a given ontology level (hypergeometric test)
///
/// ### Params
///
/// * `target_set` - The HashSet of target genes.
/// * `level` - The ontology level to process.
/// * `go_obj` The `GeneOntology` structure with the needed data.
/// * `min_genes` - Minimum number of genes in a given term.
/// * `gene_universe_length` - Size of the gene universe
/// * `elim_threshold` - Elimination threshold. Below that threshold the genes
///                      from that term are eliminated from its ancestors
/// * `debug` - Shall debug messages be printed.
///
/// ### Return
///
/// The `GoElimLevelResults` with the calculated statistics.
pub fn process_ontology_level(
    target_set: &FxHashSet<String>,
    level: &String,
    go_obj: &mut GeneOntology,
    min_genes: usize,
    gene_universe_length: usize,
    elim_threshold: f64,
) -> GoElimLevelResults {
    // Get the identifiers of that level and clean everything up
    let binding: Vec<String> = Vec::new();

    let level_ids = go_obj.get_level_ids(level).unwrap_or(&binding);
    let level_data = go_obj.get_genes_list(level_ids);

    // Filter data based on minimum gene requirement
    let level_data_final: FxHashMap<_, _> = level_data
        .iter()
        .filter(|(_, value)| value.len() >= min_genes)
        .map(|(key, value)| (key.clone(), value))
        .collect();

    let trials = target_set.len();
    let size = level_data_final.len();

    let mut go_ids = Vec::with_capacity(size);
    let mut hits_vec = Vec::with_capacity(size);
    let mut pvals = Vec::with_capacity(size);
    let mut odds_ratios = Vec::with_capacity(size);
    let mut gene_set_lengths = Vec::with_capacity(size);

    for (key, value) in level_data_final {
        let gene_set_length = value.len();
        let hits = target_set.intersection(value).count();
        let q = hits as i64 - 1;
        let pval = if q > 0 {
            hypergeom_pval(
                q as usize,
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

    for term in go_to_remove.iter() {
        if let Some(ancestors) = go_obj.get_ancestors(term) {
            let ancestors_final: Vec<String> = ancestors.to_vec();
            if let Some(genes_to_remove) = go_obj.get_genes(term) {
                let genes_to_remove = genes_to_remove.clone();

                go_obj.remove_genes(&ancestors_final, &genes_to_remove);
            }
        }
    }

    res
}

/////////////////////
// Continuous test //
/////////////////////

/// Process a given ontology level (simple fgsea)
///
/// ### Params
///
/// * `stats` - The gene level statistic.
/// * `stat_name_indices` - A hashmap linking genes to index positions in stats.
/// * `level` - The ontology level to process.
/// * `go_obj` The `GeneOntology` structure with the needed data.
/// * `go_random_perms` - The `GeneOntologyRandomPerm` with the random permutations.
/// * `gsea_params` - The GSEA parameter.
/// * `elim_threshold` - Elimination threshold. Below that threshold the genes
///                      from that term are eliminated from its ancestors
///
/// ### Return
///
/// The `GoElimLevelResultsGsea` with the calculated statistics.
#[allow(clippy::too_many_arguments)]
pub fn process_ontology_level_fgsea_simple(
    stats: &[f64],
    stat_name_indices: &FxHashMap<&String, usize>,
    level: &String,
    go_obj: &mut GeneOntology,
    go_random_perms: &GeneOntologyRandomPerm,
    gsea_params: &GseaParams,
    elim_threshold: f64,
) -> Result<GoElimLevelResultsGsea> {
    // Get the identfiers of that level and clean everything up
    let binding: Vec<String> = Vec::new();
    let go_ids = go_obj.get_level_ids(level).unwrap_or(&binding);

    // BTreeMap to make sure the order is determistic
    let mut level_data_es: BTreeMap<String, GoIntermediaryRes> = BTreeMap::new();

    for go_id in go_ids {
        if let Some(genes) = go_obj.get_genes(go_id) {
            if genes.len() >= gsea_params.min_size && genes.len() <= gsea_params.max_size {
                // Convert gene names to indices in one step
                let mut indices: Vec<i32> = genes
                    .iter()
                    .filter_map(|gene| stat_name_indices.get(gene).map(|&i| i as i32))
                    .collect();

                indices.sort();

                if !indices.is_empty() {
                    let es_res =
                        calc_gsea_stats(stats, &indices, gsea_params.gsea_param, true, false);
                    let size = indices.len();
                    level_data_es.insert(go_id.clone(), (es_res.0, size, es_res.1));
                }
            }
        }
    }

    let mut pathway_scores: Vec<f64> = Vec::with_capacity(level_data_es.len());
    let mut pathway_sizes: Vec<usize> = Vec::with_capacity(level_data_es.len());
    let mut leading_edge_indices: Vec<Vec<i32>> = Vec::with_capacity(level_data_es.len());

    for v in level_data_es.values() {
        pathway_scores.push(v.0);
        pathway_sizes.push(v.1);
        leading_edge_indices.push(v.2.clone());
    }

    let level_res: GseaResults<'_> =
        go_random_perms.get_gsea_res_simple(&pathway_scores, &pathway_sizes)?;

    let go_to_remove: Vec<&String> = level_res
        .pvals
        .iter()
        .zip(level_data_es.keys())
        .filter(|(pval, _)| *pval <= &elim_threshold) // Fixed double reference and typo
        .map(|(_, go_id)| go_id)
        .collect();

    for term in go_to_remove.iter() {
        if let Some(ancestors) = go_obj.get_ancestors(term) {
            let ancestors_final: Vec<String> = ancestors.to_vec();
            if let Some(genes_to_remove) = go_obj.get_genes(term) {
                let genes_to_remove = genes_to_remove.clone();

                go_obj.remove_genes(&ancestors_final, &genes_to_remove);
            }
        }
    }

    Ok(GoElimLevelResultsGsea {
        go_ids: level_data_es.into_keys().collect(),
        es: pathway_scores.clone(),
        nes: level_res.nes,
        size: pathway_sizes.clone(),
        pvals: level_res.pvals,
        n_more_extreme: level_res.n_more_extreme,
        ge_zero: level_res.ge_zero,
        le_zero: level_res.le_zero,
        leading_edge: leading_edge_indices,
    })
}
