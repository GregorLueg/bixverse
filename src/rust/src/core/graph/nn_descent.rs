#![allow(dead_code)]

use faer::MatRef;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use rustc_hash::FxHashSet;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::core::graph::knn::AnnDist;

/// Packed neighbours
///
/// ### Fields
///
/// * `pid` - Id
/// * `distance` - Distance
/// * `is_new` - 1 or 0 to identify if new. Using using `u32` for better memory
///   alignment
#[repr(C)]
#[derive(Clone, Copy, Debug)]
struct Neighbour {
    pid: u32,
    dist: f32,
    is_new: u32,
}

impl Neighbour {
    /// Generate a new neighbour
    ///
    /// ### Params
    ///
    /// * `pid` - Pair id
    /// * `dist` - Distance
    /// * `is_new` - Is this a new neighbour
    ///
    /// ### Returns
    ///
    /// Returns the initialised self
    #[inline(always)]
    fn new(pid: usize, dist: f32, is_new: bool) -> Self {
        Self {
            pid: pid as u32,
            dist,
            is_new: is_new as u32,
        }
    }

    /// Getter to see if the cell is new.
    ///
    /// ### Returns
    ///
    ///
    #[inline(always)]
    fn is_new(&self) -> bool {
        self.is_new != 0
    }

    /// Getter to see if the cell is new.
    ///
    /// ### Returns
    ///
    /// Returns the pair id
    #[inline(always)]
    fn pid(&self) -> usize {
        self.pid as usize
    }
}

/// NN-Descent graph builder
pub struct NNDescent {
    vectors_flat: Vec<f32>,
    dim: usize,
    n: usize,
    metric: AnnDist,
}

impl NNDescent {
    pub fn build() {}

    fn run() {}

    fn initialise_ranomd() {}

    fn local_join() {}

    fn update_neighbours() {}

    /// Fast distance calculation with unsafe pointer arithmetic
    ///
    /// ### Params
    ///
    /// * `i` - Sample index i
    /// * `j` - Sample index j
    ///
    /// ### Returns
    ///
    /// The distance between the two samples
    #[inline(always)]
    unsafe fn distance(&self, i: usize, j: usize) -> f32 {
        let ptr_i = self.vectors_flat.as_ptr().add(i * self.dim);
        let ptr_j = self.vectors_flat.as_ptr().add(j * self.dim);

        match self.metric {
            AnnDist::Euclidean => {
                let mut sum = 0_f32;
                for k in 0..self.dim {
                    let diff = *ptr_i.add(k) - *ptr_j.add(k);
                    sum += diff * diff;
                }
                // no square rooting needed
                sum
            }
            AnnDist::Cosine => {
                let mut dot = 0_f32;
                let mut norm_a = 0_f32;
                let mut norm_b = 0_f32;

                for k in 0..self.dim {
                    let a = *ptr_i.add(k);
                    let b = *ptr_j.add(k);
                    dot += a * b;
                    norm_a += a * a;
                    norm_b += b * b;
                }

                1_f32 - (dot / (norm_a.sqrt() * norm_b.sqrt()))
            }
        }
    }
}
