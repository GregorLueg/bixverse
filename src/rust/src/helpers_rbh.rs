use faer::Mat;
use std::collections::{BTreeMap, HashSet};

use crate::utils_rust::*;
use crate::utils_stats::set_similarity;

pub struct RbhTripletStruc<'a> {
    pub t1: &'a str,
    pub t2: &'a str,
    pub sim: f64,
}

/// Calculates the reciprocal best hits based on set similarities.
pub fn calculate_rbh_set<'a>(
    origin_modules: &'a BTreeMap<String, HashSet<String>>,
    target_modules: &'a BTreeMap<String, HashSet<String>>,
    overlap_coefficient: bool,
    min_similarity: f64,
    debug: bool,
) -> Vec<RbhTripletStruc<'a>> {
    let similarities: Vec<Vec<f64>> = origin_modules
        .values()
        .map(|module_i| {
            target_modules
                .values()
                .map(|module_l| {
                    let sim = set_similarity(module_i, module_l, overlap_coefficient);
                    sim
                })
                .collect()
        })
        .collect();

    let names_origin: Vec<&String> = origin_modules.keys().collect();

    let names_targets: Vec<&String> = target_modules.keys().collect();

    let mat_data: Vec<f64> = flatten_vector(similarities);

    let max_sim = array_f64_max(&mat_data);

    if debug {
        println!(
            "The observed max_similarity was {}.\nThe values are {:?}",
            max_sim, mat_data,
        );
    }

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

        let sim_mat = Mat::from_fn(nrow, ncol, |j, i| mat_data[i + j * ncol]);

        let row_maxima: Vec<f64> = sim_mat
            .row_iter()
            .map(|x| {
                let row: Vec<f64> = x.iter().cloned().collect();
                array_f64_max(&row)
            })
            .collect();

        let col_maxima: Vec<f64> = sim_mat
            .col_iter()
            .map(|x| {
                let col: Vec<f64> = x.iter().cloned().collect();
                array_f64_max(&col)
            })
            .collect();

        let mut matching_pairs: Vec<RbhTripletStruc> = Vec::new();

        for r in 0..nrow {
            for c in 0..ncol {
                let value = sim_mat[(r, c)];
                if value == row_maxima[r] && value == col_maxima[c] {
                    // matching_pairs.push((
                    //     names_origin[r].to_string(),
                    //     names_targets[c].to_string(),
                    //     value,
                    // ));
                    matching_pairs.push(RbhTripletStruc {
                        t1: names_origin[r],
                        t2: names_targets[c],
                        sim: value,
                    })
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
