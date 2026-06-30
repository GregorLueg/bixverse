use bixverse_rs::core::math::matrix_helpers::rank_matrix_col;
use bixverse_rs::enrichment::enrichment_r_wrapper::get_gsva_gs_indices;
use bixverse_rs::enrichment::singscore::*;
use bixverse_rs::prelude::*;
use extendr_api::prelude::*;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_singscore;
    fn rs_rank_matrix_col;
    fn rs_rank_matrix_col_stable;
    fn rs_singscore_single;
    fn rs_singscore_multi;
    fn rs_singscore_permutation_test;
}

/////////////
// Helpers //
/////////////

/// Parse the rank type
///
/// ### Params
///
/// * `stable` - Shall the stable ranks based on a reference set of genes be
///   used.
///
/// ### Returns
///
/// The [SingscoreRankType]
fn parse_rank_type(stable: bool) -> SingscoreRankType {
    if stable {
        SingscoreRankType::Stable
    } else {
        SingscoreRankType::Standard
    }
}

/// Parse optional indices
///
/// ### Params
///
/// * `obj` - The R memory object
///
/// ### Returns
///
/// The optional indices (shifted to 0-index)
fn parse_optional_indices(obj: Robj) -> extendr_api::Result<Option<Vec<usize>>> {
    if obj.is_null() {
        Ok(None)
    } else {
        let v: Vec<i32> = obj.try_into()?;
        Ok(Some(v.r_int_convert_shift()))
    }
}

/// Transform an optional vector to an R object helper
///
/// ### Params
///
/// * `v` -
///
/// ### Returns
///
/// The underlying R object
fn optional_vec_to_robj(v: Option<Vec<f64>>) -> Robj {
    match v {
        Some(vec) => Robj::from(vec),
        None => r!(NULL),
    }
}

///////////////
// Functions //
///////////////

/// Gene rank matrix
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Ranks the columns of a matrix with the tie method.
///
/// @param exp Numerical matrix. The expression matrix (rows = genes, columns =
/// samples).
///
/// @return Returns a matrix of ranks with the same shape as `exp`.
///
/// @export
#[extendr]
fn rs_rank_matrix_col(exp: RMatrix<f64>) -> RArray<f64, 2> {
    let exp = r_matrix_to_faer(&exp);
    let result = rank_matrix_col(&exp);

    faer_to_r_matrix(result.as_ref())
}

/// Stable-gene rank matrix
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Unit-normalised column ranks computed against a set of stable genes. Use
/// with the `stable = TRUE` setting of the singscore functions.
///
/// @param exp Numerical matrix. The expression matrix (rows = genes, columns =
/// samples).
/// @param stable_gene_indices Integer vector of stable genes. One-indexed.
///
/// @return Returns a matrix of normalised ranks with the same shape as `exp`.
///
/// @export
#[extendr]
fn rs_rank_matrix_col_stable(
    exp: RMatrix<f64>,
    stable_gene_indices: Vec<i32>,
) -> extendr_api::Result<RArray<f64, 2>> {
    let exp = r_matrix_to_faer(&exp);
    let stable: Vec<usize> = stable_gene_indices.r_int_convert_shift();
    let result = rank_matrix_col_stable(&exp, &stable);
    Ok(faer_to_r_matrix(result.as_ref()))
}

/// Rust version of singscore for a single gene set
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Rust-based implementation of singscore for a single up-regulated gene set
/// with an optional paired down-regulated set.
///
/// @param ranks Numerical matrix. The ranked expression matrix (rows = genes,
/// columns = samples). Produce with column-wise ranks (standard) or with
/// [bixverse::rs_rank_matrix_col_stable()] (stable).
/// @param up_set Integer vector. One-indexed gene indices of the up set.
/// @param down_set Integer vector or NULL. One-indexed gene indices of the
/// optional down set.
/// @param center_score Boolean. Centre scores around 0. Disabled internally
/// when `known_direction = FALSE`.
/// @param known_direction Boolean. Whether the up-set direction is known.
/// Becomes irrelevant when `down_set` is also provided.
/// @param stable Boolean. If `TRUE`, use stable-gene score bounds.
///
/// @return A named list with `TotalScore`, `TotalDispersion`, and (when
/// `down_set` is provided) `UpScore`, `UpDispersion`, `DownScore`,
/// `DownDispersion`.
///
/// @export
#[extendr]
fn rs_singscore_single(
    ranks: RMatrix<f64>,
    up_set: Vec<i32>,
    down_set: Robj,
    center_score: bool,
    known_direction: bool,
    stable: bool,
) -> extendr_api::Result<List> {
    let ranks_faer = r_matrix_to_faer(&ranks);
    let up: Vec<usize> = up_set.r_int_convert_shift();
    let down = parse_optional_indices(down_set)?;

    let result = singscore_single(
        &ranks_faer,
        &up,
        down.as_deref(),
        center_score,
        known_direction,
        parse_rank_type(stable),
    );

    let mut out = List::new(6);
    out.set_elt(0, Robj::from(result.total_score))?;
    out.set_elt(1, Robj::from(result.total_dispersion))?;
    out.set_elt(2, optional_vec_to_robj(result.up_score))?;
    out.set_elt(3, optional_vec_to_robj(result.up_dispersion))?;
    out.set_elt(4, optional_vec_to_robj(result.down_score))?;
    out.set_elt(5, optional_vec_to_robj(result.down_dispersion))?;
    out.set_names(vec![
        "total_score",
        "total_dispersion",
        "up_score",
        "up_dispersion",
        "down_score",
        "down_dispersion",
    ])?;

    Ok(out)
}

/// Rust version of singscore for many gene sets
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Rust-based implementation of singscore over many gene sets with optional
/// paired down sets.
///
/// @param ranks Numerical matrix. The ranked expression matrix.
/// @param up_list List. Up gene sets as zero-indexed indices. See
/// [bixverse::rs_prepare_gsva_gs()].
/// @param down_list List or NULL. Paired down gene sets, same length and
/// ordering as `up_list`.
/// @param center_score Boolean.
/// @param known_direction Boolean.
/// @param stable Boolean.
///
/// @return A named list with
/// \itemize{
///   \item `scores` - Numerical matrix with the scores
///   \item `dispersion` - Numerical matrix with the dispersions
/// }
/// Both matrices are of shape gene_sets × samples.
///
/// @export
#[extendr]
fn rs_singscore_multi(
    ranks: RMatrix<f64>,
    up_list: List,
    down_list: Robj,
    center_score: bool,
    known_direction: bool,
    stable: bool,
) -> extendr_api::Result<List> {
    let ranks_faer = r_matrix_to_faer(&ranks);
    let up_sets = get_gsva_gs_indices(up_list)?;

    let down_sets: Option<Vec<Vec<usize>>> = if down_list.is_null() {
        None
    } else {
        let list: List = down_list.try_into()?;
        Some(get_gsva_gs_indices(list)?)
    };

    let result = singscore_multi(
        &ranks_faer,
        &up_sets,
        down_sets.as_deref(),
        center_score,
        known_direction,
        parse_rank_type(stable),
    );

    let mut out = List::new(2);
    out.set_elt(0, Robj::from(faer_to_r_matrix(result.scores.as_ref())))?;
    out.set_elt(1, Robj::from(faer_to_r_matrix(result.dispersions.as_ref())))?;
    out.set_names(vec!["scores", "dispersions"])?;

    Ok(out)
}

/// Rust version of the singscore permutation test
///
/// @description
/// `r lifecycle::badge("experimental")`
/// For `n_permutations` iterations, draws random gene indices of the same
/// total size as the real gene set(s), scores them with the same options as
/// the real call, and builds a per-sample null distribution. Returns empirical
/// one-tailed p-values: `max(1 / n_permutations, mean(null > observed))`.
///
/// @param ranks Numerical matrix. The ranked expression matrix.
/// @param up_set Integer vector. Zero-indexed.
/// @param down_set Integer vector or NULL.
/// @param center_score,known_direction,stable Booleans. Should match the
/// values used for the real [bixverse::rs_singscore_single()] call.
/// @param n_permutations Integer. Number of random draws (B).
/// @param seed Integer. RNG seed.
///
/// @return A named list with `observed_scores` (length n_samples),
/// `null_distribution` (B × n_samples matrix), and `p_values`
/// (length n_samples).
///
/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_singscore_permutation_test(
    ranks: RMatrix<f64>,
    up_set: Vec<i32>,
    down_set: Robj,
    center_score: bool,
    known_direction: bool,
    stable: bool,
    n_permutations: usize,
    seed: usize,
) -> extendr_api::Result<List> {
    let ranks_faer = r_matrix_to_faer(&ranks);
    let up: Vec<usize> = up_set.iter().map(|&x| x as usize).collect();
    let down = parse_optional_indices(down_set)?;

    let result = singscore_permutation_test(
        &ranks_faer,
        &up,
        down.as_deref(),
        center_score,
        known_direction,
        parse_rank_type(stable),
        n_permutations,
        seed as u64,
    );

    let mut out = List::new(3);
    out.set_elt(0, Robj::from(result.observed_scores))?;
    out.set_elt(
        1,
        Robj::from(faer_to_r_matrix(result.null_distribution.as_ref())),
    )?;
    out.set_elt(2, Robj::from(result.pvals))?;
    out.set_names(vec!["observed_scores", "null_distribution", "p_values"])?;

    Ok(out)
}
