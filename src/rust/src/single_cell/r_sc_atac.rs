//! Single cell ATAC: the R bindings of `bixverse-atac`.
//!
//! Mirrors the Python bindings. Stores are written by path, regions travel as
//! parallel `chrom`, `start`, `end` vectors (0-based, half-open), cells and
//! features as 0-based store indices. The R wrappers hold the validation.

use std::path::{Path, PathBuf};

use bixverse_r::prelude::*;
use bixverse_r::single_cell::atac::errors::AtacErrors;
use bixverse_r::single_cell::atac::integration::rlsi::{rlsi_on_store, RlsiParams};
use bixverse_r::single_cell::atac::io::annotation::{
    closest_genes, gene_activity_regions, read_gtf_genes, tss_sites, GeneActivityParams, GeneRecord,
};
use bixverse_r::single_cell::atac::io::fragments::{barcode_fragment_counts, Coverage};
use bixverse_r::single_cell::atac::io::merge::{merge_fragment_files, FragmentInput};
use bixverse_r::single_cell::atac::io::tenx_peaks::tenx_peaks_to_store;
use bixverse_r::single_cell::atac::peaks::calling::{
    call_peaks, parse_tag_source, KeepDup, NarrowPeak, PeakCallParams,
};
use bixverse_r::single_cell::atac::peaks::regions::{
    genome_tiles, overlaps_any, read_bed, reduce_regions, region,
};
use bixverse_r::single_cell::atac::processing::amulet::{amulet, AmuletParams};
use bixverse_r::single_cell::atac::processing::counting::{
    count_fragments, parse_count_mode, write_counts_store,
};
use bixverse_r::single_cell::atac::processing::features::{
    feature_total_stats, parse_top_feature_cutoff, select_top_features, TopFeatureCutoff,
};
use bixverse_r::single_cell::atac::processing::links::{peak_to_gene_links, LinkParams};
use bixverse_r::single_cell::atac::processing::lsi::{
    lsi_on_store, project_lsi, read_tf_idf, set_tf_idf_on_store, LsiParams,
};
use bixverse_r::single_cell::atac::processing::profiles::fragment_profiles;
use bixverse_r::single_cell::atac::processing::qc::{fragment_qc, FragmentQcParams};
use bixverse_r::single_cell::sc_data::store::interval_index::FeatureInterval;
use bixverse_r::single_cell::sc_data::store::reader::ParallelSparseReader;
use bixverse_r::single_cell::sc_data::store::v5::parse_feature_kind;
use bixverse_r::single_cell::sc_processing::pca::{parse_pca_solver, PcaSolver};
use extendr_api::prelude::*;
use extendr_api::{Error, Result};
use faer::Mat;

use crate::single_cell::utils::max_nnz;

////////////////////
// extendr Module //
////////////////////

extendr_module! {
    mod r_sc_atac;
    fn rs_atac_read_bed;
    fn rs_atac_genome_tiles;
    fn rs_atac_reduce_regions;
    fn rs_atac_overlaps_any;
    fn rs_atac_tss_positions;
    fn rs_atac_gene_activity_regions;
    fn rs_atac_closest_genes;
    fn rs_atac_tenx_peaks_to_store;
    fn rs_atac_fragments_to_store;
    fn rs_atac_fragment_qc;
    fn rs_atac_fragment_profiles;
    fn rs_atac_barcode_fragment_counts;
    fn rs_atac_merge_fragments;
    fn rs_atac_peak_to_gene;
    fn rs_atac_call_peaks;
    fn rs_atac_amulet;
    fn rs_atac_set_tf_idf;
    fn rs_atac_feature_totals;
    fn rs_atac_lsi;
    fn rs_atac_project_lsi;
    fn rs_atac_rlsi;
}

/////////////
// Helpers //
/////////////

/// Convert an ATAC error into an extendr error
///
/// ### Params
///
/// * `e` - The error
///
/// ### Returns
///
/// The extendr error carrying the message.
fn atac_err(e: AtacErrors) -> extendr_api::Error {
    extendr_api::Error::Other(e.to_string())
}

/// The SVD solver of the LSIs by name
///
/// ### Params
///
/// * `solver` - `"auto"`, `"exact"`, `"randomised"` or `"covariance"`
///
/// ### Returns
///
/// The [PcaSolver], or an error for an unknown name.
fn lsi_solver(solver: &str) -> Result<PcaSolver> {
    parse_pca_solver(solver).ok_or_else(|| {
        Error::Other(format!(
            "solver must be 'auto', 'exact', 'randomised' or 'covariance', got '{solver}'"
        ))
    })
}

/// 0-based R integers as `usize` indices
///
/// ### Params
///
/// * `x` - The R integers, 0-based
///
/// ### Returns
///
/// The indices, or an error for a negative value.
fn index_vec(x: &[i32]) -> Result<Vec<usize>> {
    x.iter()
        .map(|&v| usize::try_from(v).map_err(|_| Error::Other(format!("negative index {v}"))))
        .collect()
}

/// 0-based R integers as `u32` values
///
/// ### Params
///
/// * `x` - The R integers
///
/// ### Returns
///
/// The values, or an error for a negative value.
fn u32_vec(x: &[i32]) -> Result<Vec<u32>> {
    x.iter()
        .map(|&v| u32::try_from(v).map_err(|_| Error::Other(format!("negative value {v}"))))
        .collect()
}

/// Regions from parallel R vectors
///
/// ### Params
///
/// * `chrom` - Chromosomes
/// * `start` - 0-based starts
/// * `end` - Exclusive ends
///
/// ### Returns
///
/// The regions, or an error for vectors of different lengths, a negative
/// coordinate or an empty region.
fn regions_from(chrom: &[String], start: &[i32], end: &[i32]) -> Result<Vec<FeatureInterval>> {
    if chrom.len() != start.len() || chrom.len() != end.len() {
        return Err(Error::Other(
            "chrom, start and end must have the same length".into(),
        ));
    }
    let start = u32_vec(start)?;
    let end = u32_vec(end)?;
    chrom
        .iter()
        .zip(start)
        .zip(end)
        .map(|((c, s), e)| {
            if e <= s {
                Err(Error::Other(format!("empty region {c}:{s}-{e}")))
            } else {
                Ok(region(c, s, e))
            }
        })
        .collect()
}

/// Regions as an R list of `chrom`, `start`, `end`
///
/// ### Params
///
/// * `regions` - The regions
///
/// ### Returns
///
/// The list, coordinates as integers.
fn regions_list(regions: &[FeatureInterval]) -> List {
    list!(
        chrom = regions.iter().map(|r| r.chrom.clone()).collect::<Vec<_>>(),
        start = regions.iter().map(|r| r.start as i32).collect::<Vec<i32>>(),
        end = regions.iter().map(|r| r.end as i32).collect::<Vec<i32>>()
    )
}

/// Paired TSS chromosomes and positions
///
/// ### Params
///
/// * `chrom` - TSS chromosomes
/// * `pos` - TSS positions, 0-based
///
/// ### Returns
///
/// The TSS as `(chrom, pos)`, or an error for vectors of different lengths.
fn tss_pairs(chrom: Vec<String>, pos: &[i32]) -> Result<Vec<(String, u32)>> {
    if chrom.len() != pos.len() {
        return Err(Error::Other(
            "tss_chrom and tss_pos differ in length".into(),
        ));
    }
    Ok(chrom.into_iter().zip(u32_vec(pos)?).collect())
}

/// Coverage of a fragment pass as an R list
///
/// ### Params
///
/// * `cov` - What the pass saw
/// * `targets` - Chromosomes that carry targets
///
/// ### Returns
///
/// List of `cells_missing`, `off_target_fraction` and `off_target_chroms`
/// (the three largest chromosomes without targets).
fn coverage_list(cov: &Coverage, targets: &[&str]) -> List {
    let mut off: Vec<(&String, &u64)> = cov
        .per_chrom
        .iter()
        .filter(|(c, _)| !targets.is_empty() && !targets.contains(&c.as_str()))
        .collect();
    off.sort_by(|a, b| b.1.cmp(a.1).then(a.0.cmp(b.0)));
    list!(
        cells_missing = cov.cells_missing as f64,
        off_target_fraction = cov.off_target_fraction(),
        off_target_chroms = off
            .iter()
            .take(3)
            .map(|(c, _)| (*c).clone())
            .collect::<Vec<String>>()
    )
}

/// Flattened rows as an R matrix
///
/// ### Params
///
/// * `rows` - The rows, all the same length
/// * `ncol` - Number of columns
///
/// ### Returns
///
/// A numeric matrix of `rows.len()` x `ncol`.
fn rows_to_r_matrix(rows: &[Vec<u64>], ncol: usize) -> RArray<f64, 2> {
    let mat = Mat::from_fn(rows.len(), ncol, |i, j| rows[i][j] as f64);
    faer_to_r_matrix(mat.as_ref())
}

/////////////
// Regions //
/////////////

/// Read a BED file
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Reads the first three columns of a BED file, plain or gzip.
///
/// @param path String. The BED file.
///
/// @returns A list with `chrom`, `start` (0-based) and `end` (exclusive).
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_atac_read_bed(path: &str) -> Result<List> {
    let regions = read_bed(Path::new(path)).map_err(atac_err)?;
    Ok(regions_list(&regions))
}

/// Fixed-width tiles over chromosomes
///
/// @description
/// `r lifecycle::badge("experimental")`
///
/// @param chroms Character vector. Chromosome names.
/// @param sizes Integer vector. Their lengths.
/// @param width Integer. The tile width.
///
/// @returns A list with `chrom`, `start` and `end`.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_atac_genome_tiles(chroms: Vec<String>, sizes: &[i32], width: i32) -> Result<List> {
    if chroms.len() != sizes.len() {
        return Err(Error::Other("chroms and sizes differ in length".into()));
    }
    let width = u32::try_from(width).map_err(|_| Error::Other("negative width".into()))?;
    let cs: Vec<(String, u32)> = chroms.into_iter().zip(u32_vec(sizes)?).collect();
    Ok(regions_list(&genome_tiles(&cs, width)))
}

/// Merge overlapping and abutting regions
///
/// @description
/// `r lifecycle::badge("experimental")`
///
/// @param chrom Character vector. Chromosomes.
/// @param start Integer vector. 0-based starts.
/// @param end Integer vector. Exclusive ends.
///
/// @returns A list with `chrom`, `start` and `end` of the merged regions and
/// `members`, a list with the input indices (0-based) each one absorbed.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_atac_reduce_regions(chrom: Vec<String>, start: &[i32], end: &[i32]) -> Result<List> {
    let regions = regions_from(&chrom, start, end)?;
    let (merged, members) = reduce_regions(&regions);
    let members: Vec<Robj> = members
        .into_iter()
        .map(|m| {
            m.into_iter()
                .map(|i| i as i32)
                .collect::<Vec<i32>>()
                .into_robj()
        })
        .collect();
    let r = regions_list(&merged);
    Ok(list!(
        chrom = r.dollar("chrom")?,
        start = r.dollar("start")?,
        end = r.dollar("end")?,
        members = List::from_values(members)
    ))
}

/// Flag regions that overlap any of another set
///
/// @description
/// `r lifecycle::badge("experimental")`
///
/// @param chrom,start,end The regions to test.
/// @param other_chrom,other_start,other_end The set to test against.
///
/// @returns Logical vector, one flag per region.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_atac_overlaps_any(
    chrom: Vec<String>,
    start: &[i32],
    end: &[i32],
    other_chrom: Vec<String>,
    other_start: &[i32],
    other_end: &[i32],
) -> Result<Vec<bool>> {
    let regions = regions_from(&chrom, start, end)?;
    let other = regions_from(&other_chrom, other_start, other_end)?;
    Ok(overlaps_any(&regions, &other))
}

////////////////
// Annotation //
////////////////

/// TSS positions from a GTF
///
/// @description
/// `r lifecycle::badge("experimental")`
///
/// @param path String. The GTF, plain or gzip.
/// @param biotypes Character vector. Biotypes to keep, empty for all.
///
/// @returns A list with `chrom`, `pos` (0-based), `strand` (`-1` minus, `1`
/// plus), `gene_id` and `gene_name`, one entry per kept gene.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_atac_tss_positions(path: &str, biotypes: Vec<String>) -> Result<List> {
    let genes = read_gtf_genes(Path::new(path)).map_err(atac_err)?;
    let tss = tss_sites(&genes, &biotypes).map_err(atac_err)?;
    let kept: Vec<&GeneRecord> = genes
        .iter()
        .filter(|g| biotypes.is_empty() || biotypes.contains(&g.biotype))
        .collect();
    Ok(list!(
        chrom = tss.iter().map(|t| t.0.clone()).collect::<Vec<_>>(),
        pos = tss.iter().map(|t| t.1 as i32).collect::<Vec<i32>>(),
        strand = tss.iter().map(|t| t.2 as i32).collect::<Vec<i32>>(),
        gene_id = kept.iter().map(|g| g.gene_id.clone()).collect::<Vec<_>>(),
        gene_name = kept.iter().map(|g| g.gene_name.clone()).collect::<Vec<_>>()
    ))
}

/// Gene body plus promoter regions from a GTF
///
/// @description
/// `r lifecycle::badge("experimental")`
/// The regions of Signac's `GeneActivity()`.
///
/// @param path String. The GTF.
/// @param upstream Integer. Bases added upstream.
/// @param downstream Integer. Bases added downstream.
/// @param max_width Integer. Genes this long or longer are dropped.
/// @param biotypes Character vector. Biotypes to keep, empty for all.
///
/// @returns A list with `chrom`, `start`, `end`, `name` (unique) and
/// `gene_id` (as in the GTF).
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_atac_gene_activity_regions(
    path: &str,
    upstream: i32,
    downstream: i32,
    max_width: i32,
    biotypes: Vec<String>,
) -> Result<List> {
    let genes = read_gtf_genes(Path::new(path)).map_err(atac_err)?;
    let params = GeneActivityParams {
        upstream: upstream as u32,
        downstream: downstream as u32,
        max_width: max_width as u32,
        biotypes,
    };
    let (regions, names, ids) = gene_activity_regions(&genes, &params).map_err(atac_err)?;
    let r = regions_list(&regions);
    Ok(list!(
        chrom = r.dollar("chrom")?,
        start = r.dollar("start")?,
        end = r.dollar("end")?,
        name = names,
        gene_id = ids
    ))
}

/// The closest gene body to every region
///
/// @description
/// `r lifecycle::badge("experimental")`
///
/// @param gtf String. The GTF.
/// @param biotypes Character vector. Biotypes to keep, empty for all.
/// @param chrom,start,end The regions.
///
/// @returns A list with `gene_name` and `gene_id` (`NA` on a chromosome
/// without genes), `distance` (0 when overlapping) and `tss_distance` (region
/// centre to TSS, positive downstream).
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_atac_closest_genes(
    gtf: &str,
    biotypes: Vec<String>,
    chrom: Vec<String>,
    start: &[i32],
    end: &[i32],
) -> Result<List> {
    let regions = regions_from(&chrom, start, end)?;
    let genes: Vec<GeneRecord> = read_gtf_genes(Path::new(gtf))
        .map_err(atac_err)?
        .into_iter()
        .filter(|g| biotypes.is_empty() || biotypes.contains(&g.biotype))
        .collect();
    if genes.is_empty() {
        return Err(Error::Other(format!(
            "no gene of the GTF has a biotype in {biotypes:?}"
        )));
    }
    let hits = closest_genes(&regions, &genes);
    let pick = |f: fn(&GeneRecord) -> &String| -> Strings {
        hits.iter()
            .map(|h| match h.gene {
                Some(i) => Rstr::from(f(&genes[i]).as_str()),
                None => Rstr::na(),
            })
            .collect()
    };
    Ok(list!(
        gene_name = pick(|g| &g.gene_name),
        gene_id = pick(|g| &g.gene_id),
        distance = hits.iter().map(|h| h.distance as f64).collect::<Vec<f64>>(),
        tss_distance = hits
            .iter()
            .map(|h| h.tss_distance as f64)
            .collect::<Vec<f64>>()
    ))
}

////////////
// Ingest //
////////////

/// The `Peaks` modality of a 10x multiome h5 into a peak store
///
/// @description
/// `r lifecycle::badge("experimental")`
///
/// @param h5_path String. The Cell Ranger ARC h5.
/// @param f_path_cells String. Cell-based file to write.
/// @param f_path_genes String. Gene-based file to write.
/// @param barcodes Character vector. Cells in store order.
/// @param max_mem_gb Optional numeric. Budget of the gene file conversion.
/// @param verbose Integer. `0L` - quiet; `1L` - normal; `2L` - detailed.
///
/// @returns A list with `name`, `chrom`, `start`, `end` per peak and
/// `lib_size`, `nnz` per cell.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_atac_tenx_peaks_to_store(
    h5_path: &str,
    f_path_cells: &str,
    f_path_genes: &str,
    barcodes: Vec<String>,
    max_mem_gb: Option<f64>,
    verbose: usize,
) -> Result<List> {
    let res = tenx_peaks_to_store(
        h5_path,
        f_path_cells,
        f_path_genes,
        &barcodes,
        max_nnz(max_mem_gb),
        verbose,
    )
    .map_err(atac_err)?;
    let r = regions_list(&res.intervals);
    Ok(list!(
        name = res.peak_names,
        chrom = r.dollar("chrom")?,
        start = r.dollar("start")?,
        end = r.dollar("end")?,
        lib_size = res.lib_size.iter().map(|&x| x as f64).collect::<Vec<f64>>(),
        nnz = res.nnz.iter().map(|&x| x as f64).collect::<Vec<f64>>()
    ))
}

/// Count fragments into features and write the store
///
/// @description
/// `r lifecycle::badge("experimental")`
///
/// @param fragments String. The fragment file.
/// @param f_path_cells String. Cell-based file to write.
/// @param f_path_genes String. Gene-based file to write.
/// @param barcodes Character vector. Cells in store order.
/// @param chrom,start,end The features.
/// @param mode String. `"insertion"`, `"paired-insertion"` or `"fragment"`.
/// @param kind String. `"peak"`, `"tile"` or `"gene_activity"`.
/// @param target_size Numeric. Target size of the initial log-CPM norm.
/// @param max_mem_gb Optional numeric. Budget of the gene file conversion.
/// @param verbose Integer. `0L` - quiet; `1L` - normal; `2L` - detailed.
///
/// @returns A list with `lib_size` and `nnz` per cell, plus `cells_missing`,
/// `off_target_fraction` and `off_target_chroms` of the pass.
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_atac_fragments_to_store(
    fragments: &str,
    f_path_cells: &str,
    f_path_genes: &str,
    barcodes: Vec<String>,
    chrom: Vec<String>,
    start: &[i32],
    end: &[i32],
    mode: &str,
    kind: &str,
    target_size: f64,
    max_mem_gb: Option<f64>,
    verbose: usize,
) -> Result<List> {
    let mode = parse_count_mode(mode).ok_or_else(|| {
        Error::Other(format!(
            "mode must be 'insertion', 'paired-insertion' or 'fragment', got '{mode}'"
        ))
    })?;
    let kind = parse_feature_kind(kind)
        .ok_or_else(|| Error::Other(format!("unknown feature kind '{kind}'")))?;
    let features = regions_from(&chrom, start, end)?;
    let counts =
        count_fragments(Path::new(fragments), &barcodes, &features, mode).map_err(atac_err)?;
    let nnz: Vec<f64> = counts.indices.iter().map(|v| v.len() as f64).collect();
    let lib = write_counts_store(
        &counts,
        &features,
        &barcodes,
        kind,
        target_size as f32,
        f_path_cells,
        f_path_genes,
        max_nnz(max_mem_gb),
        verbose,
    )
    .map_err(atac_err)?;
    let targets: Vec<&str> = chrom.iter().map(String::as_str).collect();
    let cov = coverage_list(&counts.coverage, &targets);
    Ok(list!(
        lib_size = lib.iter().map(|&x| x as f64).collect::<Vec<f64>>(),
        nnz = nnz,
        cells_missing = cov.dollar("cells_missing")?,
        off_target_fraction = cov.dollar("off_target_fraction")?,
        off_target_chroms = cov.dollar("off_target_chroms")?
    ))
}

////////
// QC //
////////

/// Per-cell QC from a fragment file
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Signac's fragment-based metrics: fragments per cell, nucleosome signal,
/// TSS enrichment and the mitochondrial fraction.
///
/// @param fragments String. The fragment file.
/// @param barcodes Character vector. Cells, output order.
/// @param tss_chrom Character vector. TSS chromosomes.
/// @param tss_pos Integer vector. TSS positions (0-based).
/// @param mito_chroms Character vector. Names of the mitochondrial genome.
///
/// @returns A list with per-cell `n_fragments`, `nucleosome_signal`,
/// `tss_enrichment` and `mito_fraction`, plus the coverage keys (off target:
/// on a chromosome without a TSS).
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_atac_fragment_qc(
    fragments: &str,
    barcodes: Vec<String>,
    tss_chrom: Vec<String>,
    tss_pos: &[i32],
    mito_chroms: Vec<String>,
) -> Result<List> {
    let tss = tss_pairs(tss_chrom, tss_pos)?;
    let params = FragmentQcParams {
        mito_chroms,
        ..Default::default()
    };
    let qc = fragment_qc(Path::new(fragments), &barcodes, &tss, &params).map_err(atac_err)?;
    let targets: Vec<&str> = tss.iter().map(|t| t.0.as_str()).collect();
    let cov = coverage_list(&qc.coverage, &targets);
    let f64s = |v: Vec<f64>| v;
    Ok(list!(
        n_fragments = qc
            .n_fragments
            .iter()
            .map(|&x| x as f64)
            .collect::<Vec<f64>>(),
        nucleosome_signal = f64s(qc.nucleosome_signal()),
        tss_enrichment = f64s(qc.tss_enrichment()),
        mito_fraction = f64s(qc.mito_fraction()),
        cells_missing = cov.dollar("cells_missing")?,
        off_target_fraction = cov.dollar("off_target_fraction")?,
        off_target_chroms = cov.dollar("off_target_chroms")?
    ))
}

/// TSS profiles and fragment length histograms per group of cells
///
/// @description
/// `r lifecycle::badge("experimental")`
/// SnapATAC2's `TSS_profile` and `frag_size_distr`.
///
/// @param fragments String. The fragment file.
/// @param barcodes Character vector. Cells.
/// @param groups Integer vector. Group of each cell, `0..n_groups`.
/// @param tss_chrom Character vector. TSS chromosomes.
/// @param tss_pos Integer vector. TSS positions (0-based).
/// @param tss_strand Integer vector. TSS strands, `-1` minus, else plus.
/// @param half_width Integer. Profile half width.
/// @param max_len Integer. Longest fragment with its own bin.
///
/// @returns A list with `tss` (groups x `2 * half_width + 1`), `lengths`
/// (groups x `max_len + 1`, bin 0 the longer fragments) and the coverage keys.
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_atac_fragment_profiles(
    fragments: &str,
    barcodes: Vec<String>,
    groups: &[i32],
    tss_chrom: Vec<String>,
    tss_pos: &[i32],
    tss_strand: &[i32],
    half_width: i32,
    max_len: i32,
) -> Result<List> {
    if tss_chrom.len() != tss_pos.len() || tss_chrom.len() != tss_strand.len() {
        return Err(Error::Other(
            "tss_chrom, tss_pos and tss_strand differ in length".into(),
        ));
    }
    let tss: Vec<(String, u32, i8)> = tss_chrom
        .into_iter()
        .zip(u32_vec(tss_pos)?)
        .zip(tss_strand)
        .map(|((c, p), &s)| (c, p, if s < 0 { -1 } else { 1 }))
        .collect();
    let groups = u32_vec(groups)?;
    let (half_width, max_len) = (half_width as u32, max_len as u32);
    let res = fragment_profiles(
        Path::new(fragments),
        &barcodes,
        &groups,
        &tss,
        half_width,
        max_len,
    )
    .map_err(atac_err)?;
    let targets: Vec<&str> = tss.iter().map(|t| t.0.as_str()).collect();
    let cov = coverage_list(&res.coverage, &targets);
    Ok(list!(
        tss = rows_to_r_matrix(&res.tss, 2 * half_width as usize + 1),
        lengths = rows_to_r_matrix(&res.lengths, max_len as usize + 1),
        cells_missing = cov.dollar("cells_missing")?,
        off_target_fraction = cov.dollar("off_target_fraction")?,
        off_target_chroms = cov.dollar("off_target_chroms")?
    ))
}

/// Fragments per barcode
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Counts every barcode in the file.
///
/// @param fragments String. The fragment file.
/// @param exclude_chroms Character vector. Chromosomes not counted.
///
/// @returns A list with `barcode` (sorted) and `n_fragments`.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_atac_barcode_fragment_counts(fragments: &str, exclude_chroms: Vec<String>) -> Result<List> {
    let counts =
        barcode_fragment_counts(Path::new(fragments), &exclude_chroms).map_err(atac_err)?;
    let (bc, n): (Vec<String>, Vec<f64>) = counts.into_iter().map(|(b, n)| (b, n as f64)).unzip();
    Ok(list!(barcode = bc, n_fragments = n))
}

/// Merge fragment files into one indexed file
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Barcodes get a prefix and suffix per input, so samples stay apart.
///
/// @param paths Character vector. The inputs.
/// @param prefixes Character vector. Barcode prefix per input.
/// @param suffixes Character vector. Barcode suffix per input.
/// @param out String. Output path; the index goes to `<out>.tbi`.
///
/// @returns A list with `n_fragments` per input and `chroms` of the output.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_atac_merge_fragments(
    paths: Vec<String>,
    prefixes: Vec<String>,
    suffixes: Vec<String>,
    out: &str,
) -> Result<List> {
    if paths.len() != prefixes.len() || paths.len() != suffixes.len() {
        return Err(Error::Other(
            "paths, prefixes and suffixes differ in length".into(),
        ));
    }
    let inputs: Vec<FragmentInput> = paths
        .into_iter()
        .zip(prefixes)
        .zip(suffixes)
        .map(|((path, prefix), suffix)| FragmentInput {
            path: PathBuf::from(path),
            prefix,
            suffix,
        })
        .collect();
    let res = merge_fragment_files(&inputs, Path::new(out)).map_err(atac_err)?;
    Ok(list!(
        n_fragments = res
            .n_fragments
            .iter()
            .map(|&x| x as f64)
            .collect::<Vec<f64>>(),
        chroms = res.chroms
    ))
}

///////////
// Links //
///////////

/// Peak to gene links over groups of cells
///
/// @description
/// `r lifecycle::badge("experimental")`
/// ArchR's `addPeak2GeneLinks()` on pseudobulked counts.
///
/// @param peak_counts Numeric matrix. Groups x peaks, summed counts.
/// @param gene_counts Numeric matrix. Groups x genes, summed counts.
/// @param chrom,start,end The peaks.
/// @param tss_chrom Character vector. TSS chromosome per gene column.
/// @param tss_pos Integer vector. TSS position per gene column.
/// @param max_distance Integer. Largest peak to TSS distance.
/// @param scale Numeric. Normalisation scale.
///
/// @returns A list with `peak`, `gene` (0-based column indices), `distance`,
/// `r`, `p` and `fdr`, one entry per tested pair.
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_atac_peak_to_gene(
    peak_counts: RMatrix<f64>,
    gene_counts: RMatrix<f64>,
    chrom: Vec<String>,
    start: &[i32],
    end: &[i32],
    tss_chrom: Vec<String>,
    tss_pos: &[i32],
    max_distance: i32,
    scale: f64,
) -> Result<List> {
    let peaks = regions_from(&chrom, start, end)?;
    let tss = tss_pairs(tss_chrom, tss_pos)?;
    let pc = r_matrix_to_faer(&peak_counts);
    let gc = r_matrix_to_faer(&gene_counts);
    let params = LinkParams {
        max_distance: max_distance as u32,
        scale,
    };
    let links = peak_to_gene_links(pc, gc, &peaks, &tss, &params).map_err(atac_err)?;
    Ok(list!(
        peak = links.iter().map(|l| l.peak as i32).collect::<Vec<i32>>(),
        gene = links.iter().map(|l| l.gene as i32).collect::<Vec<i32>>(),
        distance = links
            .iter()
            .map(|l| l.distance as f64)
            .collect::<Vec<f64>>(),
        r = links.iter().map(|l| l.r).collect::<Vec<f64>>(),
        p = links.iter().map(|l| l.p).collect::<Vec<f64>>(),
        fdr = links.iter().map(|l| l.fdr).collect::<Vec<f64>>()
    ))
}

//////////////////
// Peak calling //
//////////////////

/// Call peaks per group of cells
///
/// @description
/// `r lifecycle::badge("experimental")`
/// A port of MACS3 `callpeak` without control.
///
/// @param fragments String. The fragment file.
/// @param barcodes Character vector. Cells.
/// @param groups Integer vector. Group of each cell, `0..n_groups`.
/// @param tag_source String. `"insertion"` or `"fragment_start"`.
/// @param keep_dup Optional integer. Tags kept per position, `NULL` for all.
/// @param shift Integer. Tag shift before extension.
/// @param extsize Integer. Tag extension.
/// @param genome_size Numeric. Effective genome size.
/// @param llocal Integer. Local lambda window.
/// @param q_value Numeric. q-value cutoff.
/// @param min_length Optional integer. Minimum peak length, `NULL` for
/// `extsize`.
/// @param max_gap Optional integer. Largest merged gap, `NULL` for MACS's tag
/// size.
/// @param max_mem_gb Optional numeric. Tag budget of the first pass, `NULL`
/// for one chromosome at a time.
/// @param tmp_dir Optional string. Directory of the tag scratch file, `NULL`
/// for the system one.
///
/// @returns A list with `group` (0-based), `chrom`, `start`, `end`, `summit`
/// (0-based), `pileup`, `pscore`, `fold` and `qscore`. Groups in order, peaks
/// by chromosome name and position.
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_atac_call_peaks(
    fragments: &str,
    barcodes: Vec<String>,
    groups: &[i32],
    tag_source: &str,
    keep_dup: Option<i32>,
    shift: i32,
    extsize: i32,
    genome_size: f64,
    llocal: i32,
    q_value: f64,
    min_length: Option<i32>,
    max_gap: Option<i32>,
    max_mem_gb: Option<f64>,
    tmp_dir: Option<String>,
) -> Result<List> {
    let tag_source = parse_tag_source(tag_source)
        .ok_or_else(|| Error::Other(format!("unknown tag_source '{tag_source}'")))?;
    let params = PeakCallParams {
        tag_source,
        keep_dup: keep_dup.map_or(KeepDup::All, |k| KeepDup::AtMost(k as u32)),
        shift,
        extsize: extsize as u32,
        genome_size,
        llocal: llocal as u32,
        q_value,
        min_length: min_length.map(|v| v as u32),
        max_gap: max_gap.map(|v| v as u32),
        max_mem_gb,
        tmp_dir: tmp_dir.map(PathBuf::from),
    };
    let groups = u32_vec(groups)?;
    let peaks =
        call_peaks(Path::new(fragments), &barcodes, &groups, Some(params)).map_err(atac_err)?;
    let flat: Vec<(i32, &NarrowPeak)> = peaks
        .iter()
        .enumerate()
        .flat_map(|(g, ps)| ps.iter().map(move |p| (g as i32, p)))
        .collect();
    let ints =
        |f: fn(&NarrowPeak) -> u32| flat.iter().map(|(_, p)| f(p) as i32).collect::<Vec<i32>>();
    let dbls =
        |f: fn(&NarrowPeak) -> f32| flat.iter().map(|(_, p)| f(p) as f64).collect::<Vec<f64>>();
    Ok(list!(
        group = flat.iter().map(|&(g, _)| g).collect::<Vec<i32>>(),
        chrom = flat
            .iter()
            .map(|(_, p)| p.chrom.clone())
            .collect::<Vec<_>>(),
        start = ints(|p| p.start),
        end = ints(|p| p.end),
        summit = ints(|p| p.summit),
        pileup = dbls(|p| p.pileup),
        pscore = dbls(|p| p.pscore),
        fold = dbls(|p| p.fold),
        qscore = dbls(|p| p.qscore)
    ))
}

////////////////
// Multiplets //
////////////////

/// AMULET multiplet detection from fragment overlaps
///
/// @description
/// `r lifecycle::badge("experimental")`
///
/// @param fragments String. The fragment file.
/// @param barcodes Character vector. Cells to test.
/// @param restrict_beds Character vector. BED files of regions to drop
/// (repeats, segmental duplications, blacklists).
/// @param max_overlap Integer. Overlapping fragments a diploid cell can
/// explain.
/// @param max_fragment_len Integer. Longest valid fragment.
/// @param repeat_q Numeric. q-value under which a site is a repeat.
/// @param multiplet_q Numeric. q-value under which a cell is a multiplet.
/// @param chroms Optional character vector. Chromosomes to use, `NULL` for the
/// autosomes.
/// @param max_mem_gb Optional numeric. Read budget, `NULL` for one chromosome
/// at a time.
/// @param count_trailing_loci Boolean. Count loci no later fragment passes.
///
/// @returns A list with per-cell `n_overlaps`, `n_sites`, `p_value`,
/// `q_value` and `multiplet`.
///
/// @export
///
/// @references Thibodeau, et al., Genome Biol, 2021
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_atac_amulet(
    fragments: &str,
    barcodes: Vec<String>,
    restrict_beds: Vec<String>,
    max_overlap: i32,
    max_fragment_len: i32,
    repeat_q: f64,
    multiplet_q: f64,
    chroms: Nullable<Vec<String>>,
    max_mem_gb: Option<f64>,
    count_trailing_loci: bool,
) -> Result<List> {
    let params = AmuletParams {
        max_overlap: max_overlap as u32,
        max_fragment_len: max_fragment_len as u32,
        repeat_q,
        multiplet_q,
        chroms: chroms.into_option(),
        max_mem_gb,
        count_trailing_loci,
    };
    let mut restrict = Vec::new();
    for bed in &restrict_beds {
        restrict.extend(read_bed(Path::new(bed)).map_err(atac_err)?);
    }
    let res = amulet(Path::new(fragments), &barcodes, &restrict, Some(params)).map_err(atac_err)?;
    Ok(list!(
        n_overlaps = res
            .n_overlaps
            .iter()
            .map(|&x| x as i32)
            .collect::<Vec<i32>>(),
        n_sites = res.n_sites.iter().map(|&x| x as i32).collect::<Vec<i32>>(),
        p_value = res.p_value,
        q_value = res.q_value,
        multiplet = res.multiplet
    ))
}

////////////////////
// TF-IDF and LSI //
////////////////////

/// Switch a store to TF-IDF over a set of cells
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Header-only change of the norm; the chunks are not rewritten.
///
/// @param f_path_cells String. Cell-based file.
/// @param f_path_genes String. Gene-based file.
/// @param cells Optional integer vector. Cells for the IDF (0-based), `NULL`
/// for all.
/// @param scale Numeric. Scale factor.
///
/// @returns The IDF per feature.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_atac_set_tf_idf(
    f_path_cells: &str,
    f_path_genes: &str,
    cells: Nullable<Vec<i32>>,
    scale: f64,
) -> Result<Vec<f64>> {
    let cells = cells.into_option().map(|c| index_vec(&c)).transpose()?;
    let idf = set_tf_idf_on_store(f_path_cells, f_path_genes, cells.as_deref(), scale as f32)
        .map_err(atac_err)?;
    Ok(idf.iter().map(|&x| x as f64).collect())
}

/// Total counts per feature and the top features
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Signac's `FindTopFeatures()`.
///
/// @param f_path_genes String. Gene-based file.
/// @param cells Optional integer vector. Cells to sum over (0-based), `NULL`
/// for all.
/// @param min_cutoff Optional string. `"q<N>"` keeps features above the N-th
/// percentile, a number those with more counts, `NULL` or `"all"` every
/// feature.
///
/// @returns A list with `count`, `percentile` and `top_feature`.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_atac_feature_totals(
    f_path_genes: &str,
    cells: Nullable<Vec<i32>>,
    min_cutoff: Option<String>,
) -> Result<List> {
    let cutoff = match min_cutoff.as_deref() {
        None => TopFeatureCutoff::All,
        Some(s) => parse_top_feature_cutoff(s).ok_or_else(|| {
            Error::Other(format!(
                "min_cutoff must be 'qN' with N in [0, 100], a number or 'all', got '{s}'"
            ))
        })?,
    };
    let cells = cells.into_option().map(|c| index_vec(&c)).transpose()?;
    let reader = ParallelSparseReader::new(f_path_genes).to_extendr()?;
    let totals = feature_total_stats(&reader, cells.as_deref()).map_err(atac_err)?;
    let mut top = vec![false; totals.count.len()];
    for i in select_top_features(&totals, cutoff) {
        top[i] = true;
    }
    Ok(list!(
        count = totals.count,
        percentile = totals.percentile,
        top_feature = top
    ))
}

/// LSI of a TF-IDF store
///
/// @description
/// `r lifecycle::badge("experimental")`
///
/// @param f_path_cells String. Cell-based file.
/// @param f_path_genes String. Gene-based file.
/// @param cells Integer vector. Cells (0-based), output row order.
/// @param features Integer vector. Features to use (0-based).
/// @param n_comp Integer. Number of components.
/// @param solver String. `"auto"`, `"exact"`, `"randomised"` or
/// `"covariance"`.
/// @param scale_embeddings Boolean. Z-score each component.
/// @param scale_max Optional numeric. Clip of the z-scores.
/// @param seed Integer. Seed.
/// @param verbose Integer. `0L` - quiet; `1L` - normal; `2L` - detailed.
///
/// @returns A list with `embedding` (cells x `n_comp`), `loadings`
/// (features x `n_comp`), `singular_values` and `depth_cor`.
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_atac_lsi(
    f_path_cells: &str,
    f_path_genes: &str,
    cells: &[i32],
    features: &[i32],
    n_comp: usize,
    solver: &str,
    scale_embeddings: bool,
    scale_max: Option<f64>,
    seed: usize,
    verbose: usize,
) -> Result<List> {
    let solver = lsi_solver(solver)?;
    let cells = index_vec(cells)?;
    let features = index_vec(features)?;
    let params = LsiParams::new(
        n_comp,
        solver,
        scale_embeddings,
        scale_max.map(|v| v as f32),
    );
    let res = lsi_on_store(
        f_path_cells,
        f_path_genes,
        &cells,
        &features,
        &params,
        seed,
        verbose,
    )
    .map_err(atac_err)?;
    Ok(list!(
        embedding = faer_to_r_matrix(res.embedding.as_ref()),
        loadings = faer_to_r_matrix(res.loadings.as_ref()),
        singular_values = res
            .singular_values
            .iter()
            .map(|&x| x as f64)
            .collect::<Vec<f64>>(),
        depth_cor = res
            .depth_cor
            .iter()
            .map(|&x| x as f64)
            .collect::<Vec<f64>>()
    ))
}

/// Project cells of a TF-IDF store onto a fitted LSI
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Seurat's `ProjectSVD(mode = "lsi")`: the query's own TF-IDF times the
/// loadings over the singular values, no centring or scaling.
///
/// @param f_path_cells String. Cell-based file of the query.
/// @param f_path_genes String. Gene-based file of the query (TF-IDF norm).
/// @param cells Integer vector. Query cells (0-based), output row order.
/// @param features Integer vector. Query features (0-based), in the row order
/// of `loadings`.
/// @param loadings Numeric matrix. Features x components of the fit.
/// @param singular_values Numeric vector. Singular values of the fit.
///
/// @returns Numeric matrix, cells x components.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_atac_project_lsi(
    f_path_cells: &str,
    f_path_genes: &str,
    cells: &[i32],
    features: &[i32],
    loadings: RMatrix<f64>,
    singular_values: &[f64],
) -> Result<RArray<f64, 2>> {
    let cells = index_vec(cells)?;
    let features = index_vec(features)?;
    let loadings = r_matrix_to_faer_fp32(&loadings);
    let sv: Vec<f32> = singular_values.iter().map(|&x| x as f32).collect();
    let tf_idf =
        read_tf_idf(f_path_cells, f_path_genes, &cells, &features, None).map_err(atac_err)?;
    let res = project_lsi(&tf_idf, loadings.as_ref(), &sv).map_err(atac_err)?;
    Ok(faer_to_r_matrix(res.as_ref()))
}

/// rLSI integration of an embedding over a peak or tile store
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Signac's `IntegrateEmbeddings()` with `reduction = "rlsi"`: per-batch LSIs
/// give the anchor space, the anchors then correct `embedding`.
///
/// @param f_path_cells String. Cell-based file.
/// @param f_path_genes String. Gene-based file.
/// @param cells Integer vector. Cells (0-based), rows of `embedding`.
/// @param features Integer vector. Features of the per-batch LSIs (0-based).
/// @param batch Integer vector. 0-based batch per cell.
/// @param embedding Numeric matrix. Cells x dimensions to integrate.
/// @param dims Integer vector. 0-based LSI components the anchors use.
/// @param n_comp Integer. Components of the per-batch LSIs.
/// @param solver String. SVD solver.
/// @param k_anchor Integer. Neighbours for the anchors.
/// @param k_score Integer. Neighbours for the anchor scores.
/// @param k_weight Integer. Anchors per cell in the correction weights.
/// @param sd Numeric. Bandwidth divisor of the weights.
/// @param l2_norm Boolean. Scale and L2 normalise the anchor spaces.
/// @param knn_params List. The kNN parameters, see
/// [bixverse::params_knn_defaults()].
/// @param seed Integer. Seed.
/// @param verbose Integer. `0L` - quiet; `1L` - normal; `2L` - detailed.
///
/// @returns The integrated embedding, cells x dimensions.
///
/// @export
///
/// @references Stuart, et al., Nat Methods, 2021
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_atac_rlsi(
    f_path_cells: &str,
    f_path_genes: &str,
    cells: &[i32],
    features: &[i32],
    batch: &[i32],
    embedding: RMatrix<f64>,
    dims: &[i32],
    n_comp: usize,
    solver: &str,
    k_anchor: usize,
    k_score: usize,
    k_weight: usize,
    sd: f64,
    l2_norm: bool,
    knn_params: List,
    seed: usize,
    verbose: usize,
) -> Result<RArray<f64, 2>> {
    let solver = lsi_solver(solver)?;
    let params = RlsiParams {
        dims: index_vec(dims)?,
        k_anchor,
        k_score,
        k_weight,
        sd: sd as f32,
        l2_norm,
        lsi: LsiParams::new(n_comp, solver, true, None),
        knn_params: KnnParams::from_r_list(knn_params)?,
    };
    let cells = index_vec(cells)?;
    let features = index_vec(features)?;
    let batch = index_vec(batch)?;
    let embedding = r_matrix_to_faer_fp32(&embedding);
    let res = rlsi_on_store(
        f_path_cells,
        f_path_genes,
        &cells,
        &features,
        &batch,
        &embedding,
        &params,
        seed,
        verbose,
    )
    .map_err(atac_err)?;
    Ok(faer_to_r_matrix(res.as_ref()))
}
