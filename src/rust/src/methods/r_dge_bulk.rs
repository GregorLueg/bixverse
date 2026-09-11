//! Bindings for edgeR's quasi-likelihood chain (Chen, Lun and Smyth,
//! F1000Research, 2016), the limma-voom chain and the edgeR/limma building
//! blocks the bulk DGE class needs: filtering, normalisation factors, CPM,
//! voom and `removeBatchEffect`.
//!
//! The numerics are `edge-rs`, the two chains assembled in
//! [`bixverse_rs::methods::dge_bulk`]. Anything with a counts matrix of the
//! tested axis by samples goes through here, so both plain bulk and Milo's
//! neighbourhood counts use `rs_edger_ql`.

use bixverse_rs::methods::dge_bulk::{run_edger_ql, run_limma_dge, EdgeRQlParams, LimmaParams};
use bixverse_rs::methods::methods_r_wrapper::TestedFromR;
use bixverse_rs::prelude::*;
use edge_rs::core::expression::cpm;
use edge_rs::core::filtering::{filter_by_expr, FilterParams};
use edge_rs::core::normalisation::{calc_norm_factors, parse_norm_method};
use edge_rs::glm::test::Tested;
use edge_rs::limma::remove_batch_effect::remove_batch_effect;
use edge_rs::limma::voom::{voom, VoomParams};
use extendr_api::*;
use std::collections::HashMap;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_dge_bulk;
    fn rs_edger_ql;
    fn rs_limma_voom;
    fn rs_voom_normalise;
    fn rs_filter_by_expr;
    fn rs_calc_norm_factors;
    fn rs_cpm;
    fn rs_remove_batch_effect;
}

/////////////
// Helpers //
/////////////

/// Rebuilds an R matrix from a row-major buffer.
///
/// ### Params
///
/// * `data` - Row-major `n_rows * n_cols`
/// * `n_rows` - Number of rows
/// * `n_cols` - Number of columns
/// * `dimnames` - Optional dimnames to carry over from the input matrix
///
/// ### Returns
///
/// The column-major R matrix.
fn row_major_to_rmatrix(
    data: &[f64],
    n_rows: usize,
    n_cols: usize,
    dimnames: Option<Robj>,
) -> Result<RMatrix<f64>> {
    let mut out = RMatrix::new_matrix(n_rows, n_cols, |r, c| data[r * n_cols + c]);
    if let Some(dn) = dimnames {
        out.set_attrib("dimnames", dn)?;
    }
    Ok(out)
}

/// Checks that a design matrix has one row per sample.
///
/// ### Params
///
/// * `design_rows` - Rows of the design
/// * `n_samples` - Columns of the counts
///
/// ### Returns
///
/// `Ok(())` or an extendr error naming both numbers.
fn check_design_rows(design_rows: usize, n_samples: usize) -> Result<()> {
    if design_rows != n_samples {
        return Err(Error::Other(format!(
            "The design has {} rows against {} samples in the counts.",
            design_rows, n_samples
        )));
    }
    Ok(())
}

/// Reads a count matrix from R, integer or double, into a row-major buffer.
///
/// Count matrices in R are integer as often as not, and edgeR takes both.
///
/// ### Params
///
/// * `counts` - R matrix of features x samples, integer or double storage
///
/// ### Returns
///
/// The row-major values, the number of rows and the number of columns, or an
/// error if `counts` is not a numeric matrix.
fn counts_to_row_major(counts: &Robj) -> Result<(Vec<f64>, usize, usize)> {
    if let Ok(m) = RMatrix::<f64>::try_from(counts.clone()) {
        let (n_rows, n_cols) = (m.nrows(), m.ncols());
        return Ok((mat_to_flat_row_major(r_matrix_to_faer(&m)), n_rows, n_cols));
    }
    let m = RMatrix::<i32>::try_from(counts.clone())
        .map_err(|_| Error::Other("`counts` has to be a numeric matrix.".to_string()))?;
    let (n_rows, n_cols) = (m.nrows(), m.ncols());
    let raw: &[i32] = m.data();
    let mut out = vec![0.0; n_rows * n_cols];
    for c in 0..n_cols {
        for r in 0..n_rows {
            out[r * n_cols + c] = raw[c * n_rows + r] as f64;
        }
    }
    Ok((out, n_rows, n_cols))
}

///////////////
// Functions //
///////////////

/// Run the edgeR quasi-likelihood chain on a count matrix
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Runs `filterByExpr` -> `calcNormFactors` -> `glmQLFit` -> `glmQLFTest`,
/// implemented in Rust via the `edge-rs` crate and gated against edgeR 4.8.2.
/// The tested axis does not have to be genes: Milo's neighbourhood counts are
/// tested with the same call, with `filter = FALSE`.
///
/// @param counts Numeric matrix. Raw counts of features x samples. Must not
/// be normalised or log-transformed.
/// @param design Numeric matrix. The design matrix of samples x coefficients,
/// including the intercept. Needs at least two columns, since the null model
/// has to retain one.
/// @param edger_params Named list. The edgeR parameters, see
/// [bixverse::params_edger_ql()], plus either `coef` (0-indexed(!) design
/// columns to drop from the null model) or `contrast` (column-major weights
/// with `n_contrasts` columns).
///
/// @returns A list with the following elements
/// \itemize{
///   \item features_to_keep - Boolean. Which features survived the filters.
///   Spans the full feature axis of `counts`.
///   \item log_fc - Log2 fold changes of the tested coefficient or contrast.
///   \item log_cpm - Average log2 counts per million.
///   \item f_stat - The quasi-likelihood F statistic.
///   \item p_values - Raw p-values.
///   \item fdr - Benjamini-Hochberg adjusted p-values.
/// }
///
/// @references Chen, Lun and Smyth, F1000Research, 2016
///
/// @export
#[extendr]
fn rs_edger_ql(counts: RMatrix<f64>, design: RMatrix<f64>, edger_params: List) -> Result<List> {
    let n_features = counts.nrows();
    let n_samples = counts.ncols();
    let n_coef = design.ncols();
    check_design_rows(design.nrows(), n_samples)?;

    let params_map: HashMap<&str, Robj> = r_list_to_map(edger_params.clone())?;
    let tested = Tested::from_r_map(&params_map)?;
    let params = EdgeRQlParams::from_r_list(edger_params)?;

    // Both are column-major coming out of R and `run_edger_ql` reads them
    // row-major.
    let counts = mat_to_flat_row_major(r_matrix_to_faer(&counts));
    let design = mat_to_flat_row_major(r_matrix_to_faer(&design));

    let res = run_edger_ql(
        &counts, n_features, n_samples, &design, n_coef, &tested, &params,
    )
    .to_extendr()?;

    Ok(list!(
        features_to_keep = res.genes_to_keep,
        log_fc = res.log_fc,
        log_cpm = res.log_cpm,
        f_stat = res.f_stat,
        p_values = res.p_val,
        fdr = res.fdr
    ))
}

/// Run the limma linear model chain on a count matrix
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Runs optional `filterByExpr` -> `calcNormFactors` -> `voomLmFit` (or
/// limma-trend) -> `contrasts.fit` -> `eBayes` -> `topTable` for one
/// coefficient or contrast, implemented in Rust via the `edge-rs` crate and
/// gated against limma 3.66.0.
///
/// @param counts Numeric matrix. Raw counts of genes x samples. Must not be
/// normalised or log-transformed.
/// @param design Numeric matrix. The design matrix of samples x coefficients.
/// Must be full rank.
/// @param lib_size Numeric vector or NULL. Library size per sample. NULL uses
/// the column sums of `counts`. Pass the column sums from before gene filtering
/// to match edgeR, which keeps those on a subset `DGEList`.
/// @param limma_params Named list. The limma parameters, see
/// [bixverse::params_limma_voom()], plus either `coef` (a single 0-indexed(!)
/// design column) or `contrast` (column-major weights with `n_contrasts`
/// columns).
///
/// @returns A list with the following elements
/// \itemize{
///   \item features_to_keep - Boolean. Which genes survived the filters. Spans
///   the full gene axis of `counts`.
///   \item log_fc - Log2 fold changes of the tested coefficient or contrast.
///   \item ci_lower - Lower end of the 95% confidence interval on `log_fc`.
///   \item ci_upper - Upper end of the 95% confidence interval on `log_fc`.
///   \item ave_expr - Average log2 counts per million.
///   \item t_stat - Moderated t statistic.
///   \item p_values - Raw p-values.
///   \item fdr - Benjamini-Hochberg adjusted p-values.
///   \item b_stat - Log-odds of differential expression.
/// }
///
/// @references Law, et al., Genome Biol, 2014; Smyth, Stat Appl Genet Mol Biol,
/// 2004
///
/// @export
#[extendr]
fn rs_limma_voom(
    counts: Robj,
    design: RMatrix<f64>,
    lib_size: Nullable<Vec<f64>>,
    limma_params: List,
) -> Result<List> {
    let (counts, n_genes, n_samples) = counts_to_row_major(&counts)?;
    let lib_size: Option<Vec<f64>> = lib_size.into_option();
    let n_coef = design.ncols();
    check_design_rows(design.nrows(), n_samples)?;

    let params_map: HashMap<&str, Robj> = r_list_to_map(limma_params.clone())?;
    let tested = Tested::from_r_map(&params_map)?;
    let params = LimmaParams::from_r_list(limma_params)?;

    let design = mat_to_flat_row_major(r_matrix_to_faer(&design));

    let res = run_limma_dge(
        &counts,
        n_genes,
        n_samples,
        lib_size.as_deref(),
        &design,
        n_coef,
        &tested,
        &params,
    )
    .to_extendr()?;

    Ok(list!(
        features_to_keep = res.genes_to_keep,
        log_fc = res.log_fc,
        ci_lower = res.ci_lower,
        ci_upper = res.ci_upper,
        ave_expr = res.ave_expr,
        t_stat = res.t_stat,
        p_values = res.p_val,
        fdr = res.fdr,
        b_stat = res.b_stat
    ))
}

/// Voom-transform a count matrix
///
/// @description
/// `r lifecycle::badge("experimental")`
/// limma's `voom` on counts that are already filtered: log2-CPM against the
/// supplied (effective) library sizes, the mean-variance trend and the
/// precision weights. No filtering and no normalisation happen in here; pass
/// `lib.size * norm.factors` as `lib_size` to get voom on a normalised
/// DGEList.
///
/// @param counts Numeric matrix. Raw counts of genes x samples.
/// @param design Numeric matrix. The design matrix of samples x coefficients.
/// @param lib_size Numeric vector. The effective library size per sample.
/// @param span Numeric. Lowess span, only used if `adaptive_span = FALSE`.
/// @param adaptive_span Boolean. Derive the span from the number of genes, as
/// limma does since 3.56.
///
/// @returns A list with the following elements
/// \itemize{
///   \item e - Numeric matrix. The log2-CPM values, genes x samples. limma's
///   `E`.
///   \item weights - Numeric matrix. The precision weights, genes x samples.
///   \item trend_x - The mean-variance trend abscissae (average log2 count).
///   \item trend_y - The mean-variance trend ordinates (sqrt standard
///   deviation).
///   \item amean - Average log2-CPM per gene.
/// }
///
/// @references Law, et al., Genome Biol, 2014
///
/// @export
#[extendr]
fn rs_voom_normalise(
    counts: Robj,
    design: RMatrix<f64>,
    lib_size: Vec<f64>,
    span: f64,
    adaptive_span: bool,
) -> Result<List> {
    let dimnames = counts.get_attrib("dimnames");
    let (counts, n_genes, n_samples) = counts_to_row_major(&counts)?;
    let n_coef = design.ncols();
    check_design_rows(design.nrows(), n_samples)?;

    let design = mat_to_flat_row_major(r_matrix_to_faer(&design));

    let res = voom(
        &counts,
        n_genes,
        n_samples,
        &design,
        n_coef,
        Some(&lib_size),
        None,
        Some(VoomParams {
            span,
            adaptive_span,
            save_trend: true,
            ..Default::default()
        }),
    )
    .map_err(|e| Error::Other(e.to_string()))?;

    Ok(list!(
        e = row_major_to_rmatrix(&res.e, n_genes, n_samples, dimnames.clone())?,
        weights = row_major_to_rmatrix(&res.weights, n_genes, n_samples, dimnames)?,
        trend_x = res.trend_x,
        trend_y = res.trend_y,
        amean = res.amean
    ))
}

/// Filter lowly expressed genes
///
/// @description
/// `r lifecycle::badge("experimental")`
/// edgeR's `filterByExpr`, via the `edge-rs` crate.
///
/// @param counts Numeric matrix. Raw counts of genes x samples.
/// @param group Integer vector or NULL. Group per sample (e.g.
/// `as.integer(factor(x))`). If NULL, all samples form one group.
/// @param lib_size Numeric vector or NULL. Library size per sample. NULL uses
/// the column sums.
/// @param min_count Numeric. Minimum count in the median-sized library.
/// @param min_total_count Numeric. Minimum total count across all samples.
/// @param min_prop Numeric. Proportion of the smallest group beyond `large_n`
/// that has to express the gene.
///
/// @returns Boolean vector, one per gene. `TRUE` if the gene is kept.
///
/// @references Chen, Lun and Smyth, F1000Research, 2016
///
/// @export
#[extendr]
fn rs_filter_by_expr(
    counts: Robj,
    group: Nullable<Vec<i32>>,
    lib_size: Nullable<Vec<f64>>,
    min_count: f64,
    min_total_count: f64,
    min_prop: f64,
) -> Result<Vec<bool>> {
    let (counts, n_genes, n_samples) = counts_to_row_major(&counts)?;

    // 1-indexed factor codes from R to dense 0-indexed groups
    let group: Option<Vec<usize>> = match group {
        Nullable::NotNull(g) => Some(g.r_int_convert_shift()),
        Nullable::Null => None,
    };
    let lib_size: Option<Vec<f64>> = lib_size.into_option();

    filter_by_expr(
        &counts,
        n_genes,
        n_samples,
        lib_size.as_deref(),
        group.as_deref(),
        None,
        Some(FilterParams {
            min_count,
            min_total_count,
            min_prop,
            ..Default::default()
        }),
    )
    .map_err(|e| Error::Other(e.to_string()))
}

/// Calculate normalisation factors
///
/// @description
/// `r lifecycle::badge("experimental")`
/// edgeR's `calcNormFactors`, via the `edge-rs` crate.
///
/// @param counts Numeric matrix. Raw counts of genes x samples.
/// @param lib_size Numeric vector or NULL. Library size per sample. Pass the
/// pre-filter column sums after filtering genes, as edgeR keeps them. NULL uses
/// the column sums of `counts`.
/// @param norm_method String. One of
/// `c("TMM", "TMMwsp", "RLE", "upperquartile", "none")`.
///
/// @returns Numeric vector of normalisation factors, one per sample.
///
/// @references Robinson and Oshlack, Genome Biol, 2010
///
/// @export
#[extendr]
fn rs_calc_norm_factors(
    counts: Robj,
    lib_size: Nullable<Vec<f64>>,
    norm_method: &str,
) -> Result<Vec<f64>> {
    let (counts, n_genes, n_samples) = counts_to_row_major(&counts)?;
    let lib_size: Option<Vec<f64>> = lib_size.into_option();
    let method = parse_norm_method(norm_method)
        .ok_or_else(|| Error::Other(format!("Invalid normalisation method: {}", norm_method)))?;

    calc_norm_factors(
        &counts,
        n_genes,
        n_samples,
        lib_size.as_deref(),
        method,
        None,
        None,
    )
    .map_err(|e| Error::Other(e.to_string()))
}

/// Counts per million
///
/// @description
/// `r lifecycle::badge("experimental")`
/// edgeR's `cpm` on a plain count matrix, via the `edge-rs` crate.
///
/// @param counts Numeric matrix. Raw counts of genes x samples.
/// @param lib_size Numeric vector or NULL. Library size per sample, e.g.
/// `lib.size * norm.factors`. NULL uses the column sums.
/// @param log Boolean. Return log2-CPM.
/// @param prior_count Numeric. Prior count added before the log. Ignored if
/// `log = FALSE`.
///
/// @returns Numeric matrix of (log2-)CPM values, genes x samples.
///
/// @export
#[extendr]
fn rs_cpm(
    counts: Robj,
    lib_size: Nullable<Vec<f64>>,
    log: bool,
    prior_count: f64,
) -> Result<RMatrix<f64>> {
    let dimnames = counts.get_attrib("dimnames");
    let (counts, n_genes, n_samples) = counts_to_row_major(&counts)?;
    let lib_size: Option<Vec<f64>> = lib_size.into_option();

    let res = cpm(
        &counts,
        n_genes,
        n_samples,
        lib_size.as_deref(),
        None,
        log,
        prior_count,
    )
    .map_err(|e| Error::Other(e.to_string()))?;

    row_major_to_rmatrix(&res, n_genes, n_samples, dimnames)
}

/// Remove batch effects from a log-expression matrix
///
/// @description
/// `r lifecycle::badge("experimental")`
/// limma's `removeBatchEffect`, via the `edge-rs` crate. Batch gets
/// sum-to-zero contrasts, is fitted jointly with the design of interest and
/// only the batch part is subtracted. Meant for plotting and unsupervised
/// work; for testing put batch into the design.
///
/// @param x Numeric matrix. Log-expression values of genes x samples.
/// @param batch Integer vector. Batch per sample.
/// @param design Numeric matrix or NULL. The design of interest, samples x
/// coefficients, whose effects are protected. NULL is an intercept only.
///
/// @returns Numeric matrix of corrected values, genes x samples.
///
/// @references Smyth, Stat Appl Genet Mol Biol, 2004
///
/// @export
#[extendr]
fn rs_remove_batch_effect(
    x: RMatrix<f64>,
    batch: Vec<i32>,
    design: Nullable<RMatrix<f64>>,
) -> Result<RMatrix<f64>> {
    let n_genes = x.nrows();
    let n_samples = x.ncols();
    let dimnames = x.get_attrib("dimnames");
    let x = mat_to_flat_row_major(r_matrix_to_faer(&x));
    let batch: Vec<usize> = batch.r_int_convert();

    let design: Option<(Vec<f64>, usize)> = match design {
        Nullable::NotNull(d) => {
            check_design_rows(d.nrows(), n_samples)?;
            let n_coef = d.ncols();
            Some((mat_to_flat_row_major(r_matrix_to_faer(&d)), n_coef))
        }
        Nullable::Null => None,
    };

    let res = remove_batch_effect(
        &x,
        n_genes,
        n_samples,
        Some(&batch),
        None,
        None,
        design.as_ref().map(|(d, n)| (d.as_slice(), *n)),
        None,
    )
    .map_err(|e| Error::Other(e.to_string()))?;

    row_major_to_rmatrix(&res, n_genes, n_samples, dimnames)
}
