use faer::Mat;
use rustc_hash::FxHashSet;
use std::collections::BTreeMap;

use crate::utils::general::*;
use crate::utils::utils_stats::set_similarity;

/// Structure for an Rbh triplet Result
#[derive(Clone, Debug)]
pub struct RbhTripletStruc<'a> {
    pub t1: &'a str,
    pub t2: &'a str,
    pub sim: f64,
}

////////////////////
// Set similarity //
////////////////////

/// Calculates the reciprocal best hits based on set similarities.
pub fn calculate_rbh_set<'a>(
    origin_modules: &'a BTreeMap<String, FxHashSet<String>>,
    target_modules: &'a BTreeMap<String, FxHashSet<String>>,
    overlap_coefficient: bool,
    min_similarity: f64,
    debug: bool,
) -> Vec<RbhTripletStruc<'a>> {
    let names_targets: Vec<&String> = target_modules.keys().collect();
    let names_origin: Vec<&String> = origin_modules.keys().collect();

    if debug {
        println!("Target names: {:?}", names_targets)
    }

    if debug {
        println!("Origin names: {:?}", names_origin)
    }

    let similarities_flat: Vec<Vec<f64>> = origin_modules
        .values()
        .map(|v1| {
            let v1_refs: FxHashSet<&String> = v1.iter().collect();
            target_modules
                .values()
                .map(|v2| {
                    let v2_refs: FxHashSet<&String> = v2.iter().collect();
                    set_similarity(&v1_refs, &v2_refs, overlap_coefficient)
                })
                .collect()
        })
        .collect();

    let mat_data: Vec<f64> = flatten_vector(similarities_flat);

    if debug {
        println!("Flat data {:?}", mat_data)
    };

    let max_sim = array_max(&mat_data);

    let result = if max_sim < min_similarity {
        if debug {
            println!("No similarity passed the threshold.\n\n")
        }
        vec![RbhTripletStruc {
            t1: "NA",
            t2: "NA",
            sim: 0.0,
        }]
    } else {
        let nrow = names_origin.len();
        let ncol = names_targets.len();

        let sim_mat = Mat::from_fn(nrow, ncol, |i, j| mat_data[j + i * ncol]);

        if debug {
            println!("The matrix looks like: {:?}", sim_mat)
        };

        let row_maxima: Vec<&f64> = sim_mat
            .row_iter()
            .map(|x| {
                let row: Vec<&f64> = x.iter().collect();
                array_max(&row)
            })
            .collect();

        let col_maxima: Vec<&f64> = sim_mat
            .col_iter()
            .map(|x| {
                let col: Vec<&f64> = x.iter().collect();
                array_max(&col)
            })
            .collect();

        let mut matching_pairs: Vec<RbhTripletStruc> = Vec::new();

        for r in 0..nrow {
            for c in 0..ncol {
                let value = sim_mat[(r, c)];
                if value == *row_maxima[r] && value == *col_maxima[c] {
                    let triplet = RbhTripletStruc {
                        t1: names_origin[r],
                        t2: names_targets[c],
                        sim: value,
                    };

                    if debug {
                        println!("What are the matching pairs?: {:?}", triplet)
                    };

                    matching_pairs.push(triplet)
                }
            }
        }

        if debug {
            println!(
                "A total of {} RBH pairs were identified.\n\n",
                matching_pairs.len()
            );
        }

        if !matching_pairs.is_empty() {
            matching_pairs
        } else {
            vec![RbhTripletStruc {
                t1: "NA",
                t2: "NA",
                sim: 0.0,
            }]
        }
    };

    result
}
