use rand::prelude::*;
use rand::seq::SliceRandom;

pub fn split_vector_randomly(vec: Vec<f64>, x: usize, seed: u64) -> (Vec<f64>, Vec<f64>) {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut shuffled = vec.clone();
    shuffled.shuffle(&mut rng);

    let first_set = shuffled[..x].to_vec();
    let second_set = shuffled[x..].to_vec();

    (first_set, second_set)
}