use faer::{Mat, MatRef};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::assert_symmetric_mat;

///////////////////////
// KNN and Laplacian //
///////////////////////

/// Helper struct for KNN with heap
///
/// ### Fields
///
/// * `index` - Index position of that neighbours
/// * `similiarity` - Similarity value for that Neighbour
#[derive(Debug)]
struct SimilarityItem {
    index: usize,
    similarity: f64,
}

/// Equality Trait for SimilarityItem
impl Eq for SimilarityItem {}

/// `PartialEq` trait for `SimilarityItem`
///
/// Check equality between two `SimilarityItem` items in terms of similarity
///
/// ### Returns
///
/// `true` if they are the same.
impl PartialEq for SimilarityItem {
    fn eq(&self, other: &Self) -> bool {
        self.similarity == other.similarity
    }
}

/// Ord trait `SimilarityItem`
///
/// How to order `SimilarityItem`
impl Ord for SimilarityItem {
    fn cmp(&self, other: &Self) -> Ordering {
        // Reverse for min-heap (we want to keep highest similarities)
        other
            .similarity
            .partial_cmp(&self.similarity)
            .unwrap_or(Ordering::Equal)
    }
}

/// PartialOrd trait `SimilarityItem`
///
/// How to order `SimilarityItem`
impl PartialOrd for SimilarityItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Generate a KNN graph adjacency matrix from a similarity matrix
///
/// ### Params
///
/// * `similarities` - The symmetric similarity matrix.
/// * `k` - Number of neighbours to take
///
/// ### Returns
///
/// The KNN adjacency matrix
pub fn get_knn_graph_adj(similarities: &MatRef<f64>, k: usize) -> Mat<f64> {
    assert_symmetric_mat!(similarities);

    let n = similarities.nrows();
    let mut adjacency: Mat<f64> = Mat::zeros(n, n);

    // Parallelize across rows
    let rows: Vec<Vec<(usize, f64)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut heap = BinaryHeap::with_capacity(k + 1);

            // Use min-heap to keep top-k similarities
            for j in 0..n {
                if i != j {
                    let sim = similarities[(i, j)];
                    heap.push(SimilarityItem {
                        index: j,
                        similarity: sim,
                    });

                    if heap.len() > k {
                        heap.pop(); // Remove smallest
                    }
                }
            }

            heap.into_iter()
                .map(|item| (item.index, item.similarity))
                .collect()
        })
        .collect();

    // Fill adjacency matrix
    for (i, neighbors) in rows.iter().enumerate() {
        for &(j, sim) in neighbors {
            adjacency[(i, j)] = sim;
        }
    }

    // Symmetrize in parallel
    for i in 0..n {
        for j in i + 1..n {
            let val = (adjacency[(i, j)] + adjacency[(j, i)]) / 2.0;
            adjacency[(i, j)] = val;
            adjacency[(j, i)] = val;
        }
    }

    adjacency
}

/// Generate a Laplacian matrix from an adjacency matrix
///
/// ### Params
///
/// * `adjacency` - The symmetric adjacency matrix.
///
/// ### Returns
///
/// The Laplacian matrix
pub fn adjacency_to_laplacian(adjacency: &MatRef<f64>) -> Mat<f64> {
    assert_symmetric_mat!(adjacency);

    let n = adjacency.nrows();

    let mut laplacian = adjacency.cloned();

    for i in 0..n {
        let degree = adjacency.row(i).iter().sum::<f64>();
        laplacian[(i, i)] = degree - adjacency[(i, i)];

        for j in 0..n {
            if i != j {
                laplacian[(i, j)] = -adjacency[(i, j)];
            }
        }
    }

    laplacian
}

//////////////////
// Annoy search //
//////////////////

/// Tree node representation for binary space partitioning
///
/// Each node is either a split (with hyperplane and child pointers) or leaf
/// (with vector indices)
#[derive(Clone)]
enum AnnoyNode {
    /// Internal node that splits space using a hyperplane
    Split {
        /// Random hyperplane normal vector for splitting
        hyperplane: Vec<f32>,
        /// Threshold for hyperplane decision (median of projections)
        threshold: f32,
        /// Index of left child node
        left: usize,
        /// Index of right child node
        right: usize,
    },
    /// Terminal node containing actual vector indices
    Leaf {
        /// Indices of vectors that ended up in this partition
        items: Vec<usize>,
    },
}

/// Approximate Nearest Neighbors index using random binary trees
///
/// Partitions vector space using multiple random hyperplanes to enable fast
/// approximate neighbour search. Trade-off between accuracy and speed
/// controlled by number of trees. (Reminds me of Random Forests...)
///
/// ### Fields
///
/// * `trees` - Collection of binary trees, each with different random
///   partitioning.
/// * `vectors` - Original vector data for distance calculations.
/// * `n_trees` - Number of trees built (more trees = better accuracy, slower
///   queries)
pub struct AnnoyIndex {
    trees: Vec<Vec<AnnoyNode>>,
    vectors: Vec<Vec<f32>>,
    n_trees: usize,
}

impl AnnoyIndex {
    /// Creates a new Annoy index from an embedding-type matrix
    ///
    /// ### Params
    ///
    /// * `mat` - Matrix with rows = samples and columns = features.
    /// * `n_trees` - Number of random trees to build (50-100 recommended).
    /// * `seed` - Random seed for reproducible results.
    ///
    /// ### Returns
    ///
    /// Initialised AnnoyIndex ready for querying
    ///
    /// ### Algorithm Details
    ///
    /// **Tree Building Phase:**
    ///
    /// 1. For each tree, generate a random hyperplane (vector of uniform random
    ///    values in [-1,1]).
    /// 2. Project all vectors onto the hyperplane using dot products.
    /// 3. Sort projections and split at median: vectors ≤ median go left,
    ///    others go right.
    /// 4. Recursively apply steps 1-3 to each partition until ≤10 vectors
    ///    remain (leaf nodes).
    /// 5. Repeat process for n_trees independent trees using different random
    ///    seeds
    ///
    /// **Query Phase:**
    ///
    /// 1. For each tree, traverse from root following hyperplane decisions: if
    ///    query·hyperplane ≤ threshold go left, else right.
    /// 2. At boundary decisions (|query·hyperplane - threshold| < 0.5), explore
    ///    both subtrees if search budget allows.
    /// 3. Collect all vector indices from reached leaf nodes across all trees
    /// 4. Remove duplicates and compute exact Euclidean distances to all
    ///    candidates
    /// 5. Return k vectors with smallest distances
    pub fn new(mat: MatRef<f32>, n_trees: usize, seed: usize) -> Self {
        let mut rng = StdRng::seed_from_u64(seed as u64);

        let vectors: Vec<Vec<f32>> = (0..mat.nrows())
            .map(|i| mat.row(i).iter().cloned().collect())
            .collect();

        let seeds: Vec<u64> = (0..n_trees).map(|_| rng.random()).collect();

        let trees: Vec<Vec<AnnoyNode>> = seeds
            .into_par_iter()
            .map(|tree_seed| {
                let mut tree_rng = StdRng::seed_from_u64(tree_seed);
                Self::build_tree(&vectors, (0..vectors.len()).collect(), &mut tree_rng)
            })
            .collect();

        AnnoyIndex {
            trees,
            vectors,
            n_trees,
        }
    }

    /// Builds a single random tree from the given items
    ///
    /// ### Params
    ///
    /// * `vectors` - All vectors in the dataset.
    /// * `items` - Subset of vector indices to include in this tree.
    /// * `rng` - Random number generator for this tree.
    ///
    /// ### Returns
    ///
    /// Vector of AnnoyNodes representing the tree structure (index 0 = root)
    fn build_tree(vectors: &[Vec<f32>], items: Vec<usize>, rng: &mut StdRng) -> Vec<AnnoyNode> {
        let mut nodes = Vec::with_capacity(items.len() * 2);
        Self::build_node(vectors, items, &mut nodes, rng);
        nodes
    }

    /// Recursively builds tree nodes using random hyperplanes
    ///
    /// ### Params
    ///  
    /// * `vectors` - All vectors in the dataset
    /// * `items` - Items to split at this node
    /// * `nodes` - Growing list of tree nodes
    /// * `rng` - Random number generator
    ///
    /// ### Returns
    ///
    /// Index of the created node in the nodes vector.
    ///
    /// ### Implementation Details
    ///
    /// - Creates leaf if ≤10 items remain (prevents over-partitioning)
    /// - Hyperplane: random vector with components element [-1,1]  
    /// - Threshold: median of dot products (balanced split)
    /// - Fallback: creates leaf if split fails (empty partitions)
    fn build_node(
        vectors: &[Vec<f32>],
        items: Vec<usize>,
        nodes: &mut Vec<AnnoyNode>,
        rng: &mut StdRng,
    ) -> usize {
        // create leaf if too few items remain
        if items.len() <= 10 {
            let node_idx = nodes.len();
            nodes.push(AnnoyNode::Leaf { items });
            return node_idx;
        }

        let dim = vectors[0].len();
        let hyperplane: Vec<f32> = (0..dim).map(|_| rng.random_range(-1.0..1.0)).collect();

        // calculate dot products and sort to find splitting threshold
        let mut item_dots: Vec<(usize, f32)> = items
            .iter()
            .map(|&item| {
                let dot = vectors[item]
                    .iter()
                    .zip(&hyperplane)
                    .fold(0.0, |acc, (v, h)| acc + v * h);
                (item, dot)
            })
            .collect();

        item_dots.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        let threshold = item_dots[item_dots.len() / 2].1;

        // split items based on threshold
        let mut left_items = Vec::new();
        let mut right_items = Vec::new();

        for (item, dot) in item_dots {
            if dot <= threshold {
                left_items.push(item);
            } else {
                right_items.push(item);
            }
        }

        // fallback to leaf if split failed
        if left_items.is_empty() || right_items.is_empty() {
            let node_idx = nodes.len();
            nodes.push(AnnoyNode::Leaf { items });
            return node_idx;
        }

        // create split node and recursively build children
        let node_idx = nodes.len();
        nodes.push(AnnoyNode::Split {
            hyperplane: hyperplane.clone(),
            threshold,
            left: 0,
            right: 0,
        });

        let left_idx = Self::build_node(vectors, left_items, nodes, rng);
        let right_idx = Self::build_node(vectors, right_items, nodes, rng);

        // update child pointers
        if let AnnoyNode::Split {
            ref mut left,
            ref mut right,
            ..
        } = nodes[node_idx]
        {
            *left = left_idx;
            *right = right_idx;
        }

        node_idx
    }

    /// Queries the index for k nearest neighbors with enhanced search
    ///
    /// ### Params
    ///
    /// * `query_vec` - Vector to find neighbors for.
    /// * `k` - Number of neighbors to return
    /// * `search_k` - Search budget (None = k * n_trees, higher = better
    ///   recall)
    ///
    /// ### Returns
    ///
    /// Vector of k nearest neighbor indices, sorted by distance
    pub fn query(&self, query_vec: &[f32], k: usize, search_k: Option<usize>) -> Vec<usize> {
        let search_k = search_k.unwrap_or(k * self.n_trees);
        let mut candidates = Vec::with_capacity(search_k * 2);

        // query each tree with portion of search budget
        let budget_per_tree = (search_k / self.n_trees).max(k);

        for tree in &self.trees {
            self.query_tree_enhanced(tree, query_vec, &mut candidates, budget_per_tree);
        }

        // remove duplicates and calculate distances
        candidates.sort_unstable();
        candidates.dedup();

        let mut scored = Vec::with_capacity(candidates.len());
        for idx in candidates {
            let dist = self.vectors[idx]
                .iter()
                .zip(query_vec)
                .fold(0.0, |acc, (a, b)| acc + (a - b).powi(2));
            scored.push((idx, dist));
        }

        // return k closest neighbors
        scored.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        scored.into_iter().take(k).map(|(idx, _)| idx).collect()
    }

    /// Enhanced tree query with multi-path traversal and search budget
    ///
    /// ### Params
    ///
    /// * `tree` - Array of tree nodes (index 0 = root)  
    /// * `query_vec` - Query vector for similarity matching
    /// * `candidates` - Accumulating list of potential neighbor indices
    /// * `budget` - Maximum candidates to collect from this tree
    fn query_tree_enhanced(
        &self,
        tree: &[AnnoyNode],
        query_vec: &[f32],
        candidates: &mut Vec<usize>,
        budget: usize,
    ) {
        let mut remaining_budget = budget;
        Self::traverse_node_enhanced(tree, 0, query_vec, candidates, &mut remaining_budget);
    }

    /// Multi-path tree traversal with budget-controlled exploration
    ///
    /// Explores both sides of splits when query point is near the boundary,
    /// improving recall at the cost of some performance.
    ///
    /// ### Params
    ///
    /// * `tree` - Tree nodes
    /// * `node_idx` - Current node index
    /// * `query_vec` - Query vector
    /// * `candidates` - Growing candidate list
    /// * `remaining_budget` - Candidates left to collect
    fn traverse_node_enhanced(
        tree: &[AnnoyNode],
        node_idx: usize,
        query_vec: &[f32],
        candidates: &mut Vec<usize>,
        remaining_budget: &mut usize,
    ) {
        if *remaining_budget == 0 {
            return;
        }

        match &tree[node_idx] {
            AnnoyNode::Leaf { items } => {
                // Add leaf items to candidates
                let items_to_add = items.len().min(*remaining_budget);
                candidates.extend(&items[..items_to_add]);
                *remaining_budget = remaining_budget.saturating_sub(items_to_add);
            }
            AnnoyNode::Split {
                hyperplane,
                threshold,
                left,
                right,
            } => {
                let dot = query_vec
                    .iter()
                    .zip(hyperplane)
                    .fold(0.0, |acc, (q, h)| acc + q * h);

                let margin = (dot - threshold).abs();

                // Multi-path exploration when near boundary and budget allows
                if margin < 0.5 && *remaining_budget > 20 {
                    let half_budget = *remaining_budget / 2;

                    // Explore closer side first
                    let (first, second) = if dot <= *threshold {
                        (*left, *right)
                    } else {
                        (*right, *left)
                    };

                    let mut first_budget = half_budget + (*remaining_budget % 2);
                    Self::traverse_node_enhanced(
                        tree,
                        first,
                        query_vec,
                        candidates,
                        &mut first_budget,
                    );

                    let mut second_budget = half_budget;
                    Self::traverse_node_enhanced(
                        tree,
                        second,
                        query_vec,
                        candidates,
                        &mut second_budget,
                    );

                    *remaining_budget = 0;
                } else {
                    // Single-path traversal
                    let next = if dot <= *threshold { *left } else { *right };
                    Self::traverse_node_enhanced(
                        tree,
                        next,
                        query_vec,
                        candidates,
                        remaining_budget,
                    );
                }
            }
        }
    }
}
