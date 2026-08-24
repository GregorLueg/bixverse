# bixverse API index

Every documented entry point, grouped the way the package website groups them. Generated from `_pkgdown.yml` and `man/*.Rd` by `data-raw/generate_api_index.R`. Do not edit by hand.

Grep this file to check whether a function exists before calling it.

## General getters

All types of general getters that work across various classes related to co-expression module detection, graph-based clustering, etc.

- `get_metadata`: Return the metadata
- `get_params`: Get the parameters that were used.
- `get_results`: Get the final results from the class

## Gene set enrichment helpers

Everything and anything you need to do various types of gene set enrichments; hypergeometric tests, GSVA, (ss)GSEA.

- `gse_hypergeometric`: Gene set enrichment (GSE) based on a hypergeometric test.
- `gse_hypergeometric_list`: Gene set enrichment (GSE) based on a hypergeometric test over a list.
- `calc_fgsea`: Bixverse implementation of the fgsea algorithm
- `calc_fgsea_simple`: Bixverse implementation of the simple fgsea algorithm
- `calc_gsea_traditional`: Bixverse implementation of the traditional GSEA algorithm
- `calc_mitch`: Calculate a mitch gene set enrichments on contrast
- `calc_gsva`: Bixverse implementation of GSVA
- `calc_ssgsea`: Bixverse implementation of ssGSEA
- `calc_singscore`: Bixverse implementation of singscore (single gene set)
- `calc_singscore_multi`: Bixverse implementation of singscore (multiple gene sets)
- `calc_singscore_rank`: Rank an expression matrix for singscore
- `params_gsea`: Wrapper function to generate GSEA parameters
- `params_gsva`: Wrapper function to generate GSVA parameters
- `params_ssgsea`: Wrapper function to generate ssGSEA parameters

## Gene Ontology-related gene set enrichment helpers

GSE methods specifically designed with the DAG of the Gene Ontology in mind, identifying the most relevant Gene Ontology terms.

- `load_go_human_data`: Get the Gene Ontology data human
- `get_go_data_human`: Wrapper function to load and process the gene ontology data.
- `process_go_data`: Process Gene Ontology data into the right format
- `GeneOntologyElim`: Gene Ontology data
- `gse_go_elim_method`: Run gene ontology enrichment with elimination method.
- `gse_go_elim_method_list`: Run gene ontology enrichment with elimination method over a list.
- `fgsea_go_elim`: Run GO enrichment with elimination method over a continuous vectors
- `fgsea_simple_go_elim`: Run GO enrichment with elimination with fgsea simple
- `simplify_hypergeom_res`: Simplify gene set results via ontologies

## CisTarget-related helpers

Functions related to CisTarget enrichment.

- `download_cistarget_hg38`: Download CisTarget reference files for human (hg38)
- `read_motif_annotation_file`: Read in the motif annotation file
- `read_motif_ranking`: Read in the motif rankings and transform them into a matrix
- `run_cistarget`: Main function to run CisTarget
- `params_cistarget`: Wrapper function to CisTarget parameters

## Co-expression methods for bulk RNAseq

Methods to identify co-expression modules via correlations or matrix factorisations.

- `BulkCoExp`: Bulk RNAseq co-expression modules
- `preprocess_bulk_coexp`: Process the raw data
- `get_diagnostics`: Get diagnostics from a BulkModuleResult
- `get_factors`: Get factor matrices from a BulkModuleResult
- `get_nmf_gene_loadings`: Get the NMF gene loadings
- `get_nmf_modules`: Get the NMF module membership data.table
- `get_nmf_sample_activity`: Get the NMF sample activity
- `get_nmf_stability`: Get the stabilised NMF diagnostics
- `get_c_pca_factors`: Get the contrastive PCA factors
- `get_c_pca_loadings`: Get the contrastive PCA loadings
- `get_ica_stability_res`: Get the ICA component data (stability, convergence, nMI)
- `get_grid_search_res`: Get the grid search results
- `get_cor_graph`: Get correlation-based graph
- `get_diffcor_graph`: Get differential correlation-based graph
- `get_epsilon_res`: Return the epsilon data
- `get_resolution_res`: Return the resolution results
- `get_outputs`: Return the outputs
- `get_modules`: Get the module membership from a BulkModuleResult
- `cor_module_check_epsilon`: Iterate through different epsilon parameters
- `cor_module_coremo_clustering`: Generates CoReMo-based gene modules
- `cor_module_coremo_cor_sign`: Split CoReMo modules by correlation sign
- `cor_module_coremo_eigengene`: Calculate Eigengenes for CoReMo modules
- `cor_module_coremo_stability`: Assesses CoReMo-based gene module stability
- `cor_module_graph_check_res`: Iterate through Leiden resolutions for graph-based community detection.
- `cor_module_graph_final_modules`: Identify correlation-based gene modules via graphs
- `cor_module_processing`: Prepare correlation-based module detection
- `cor_module_tom`: Update the correlation matrix to a TOM
- `diffcor_module_processing`: Prepare differential correlation-based module detection
- `contrastive_pca_processing`: Prepare class for contrastive PCA
- `c_pca_plot_alphas`: Plot various alphas for the contrastive PCA
- `contrastive_pca`: Apply contrastive PCA.
- `dgrdl_grid_search`: Grid search over DGRDL parameters
- `dgrdl_result`: Run DGRDL with the specified parameters
- `ica_evaluate_comp`: Iterate over different ncomp parameters for ICA
- `ica_optimal_ncomp`: Identify stability inflection point
- `ica_processing`: Prepare class for ICA
- `ica_stabilised_results`: Run stabilised ICA with a given number of components
- `nmf_bulk`: Run non-negative matrix factorisation on a BulkCoExp
- `stabilised_nmf_bulk`: Run stabilised (multi-restart) NMF on a BulkCoExp
- `params_cor_graph`: Wrapper function for graph generation
- `params_coremo`: Wrapper function to generate CoReMo parameters
- `params_dgrdl`: Wrapper function to generate DGRDL parameters
- `params_module_membership`: Wrapper function to generate module membership parameters
- `params_ica_general`: Wrapper function for standard ICA parameters
- `params_ica_ncomp`: Wrapper function for ICA ncomp iterations
- `params_ica_randomisation`: Wrapper function for ICA randomisation

## Helpers for differential gene expression

Methods to help out with differential gene expression analyses in a structured way. Useful when you have to analyse 10's to 100's of differential gene expression results.

- `BulkDge`: Bulk RNAseq differential gene expression class
- `add_new_metadata`: Replace the meta data
- `change_gene_identifier`: Change the primary gene identifier of BulkDge
- `update_metadata_values`: Replace values in a metadata column
- `fix_meta_data_column`: Helper to fix meta-data columns to be R conform
- `remove_samples`: Remove samples from object
- `qc_bulk_dge`: QC on the bulk dge data
- `preprocess_bulk_dge`: QC on the bulk dge data (DEPRECATED!)
- `normalise_bulk_dge`: Normalise the count data for DGE.
- `batch_correction_bulk_dge`: Run a linear batch correction
- `bulk_dge_from_h5ad`: Wrapper function to generate BulkDge object from h5ad
- `calculate_dge_hedges`: Calculates the Hedge's G effect size
- `calculate_all_dges`: Calculate all possible DGE variants (DEPRECATED!)
- `calculate_dge_limma`: Calculates the Limma Voom DGE
- `calculate_pca_bulk_dge`: Calculate PCA on the expression.
- `calculate_rpkm`: RPKM calculation
- `calculate_tpm`: TPM calculation
- `run_limma_voom`: Wrapper for a Limma Voom analysis
- `hedges_g_dge`: Calculate the effect size
- `get_dge_effect_sizes`: Return the effect size results
- `get_dge_limma_voom`: Return the Limma Voom results
- `get_dge_list`: Return the DGEList
- `get_dge_qc_plot`: Return QC plots
- `get_fpkm_counts`: Return the FPKM-normalised counts
- `get_gene_lengths`: Get the gene lengths
- `get_model_fit`: Get the fitted model
- `get_tpm_counts`: Return the TPM-normalised counts

## Biomedical ontologies

For dealing with ontologies and calculating (semantic) similarities in disease, phenotype or gene ontologies.

- `OntologySim`: OntologySim class
- `pre_process_sim_onto`: Pre-process data for subsequent ontology similarity
- `calculate_information_content`: Calculate the information content for each ontology term
- `calculate_semantic_sim`: Calculate the Resnik or Lin semantic similarity
- `calculate_semantic_sim_mat`: Calculate the Resnik or Lin semantic similarity matrix
- `calculate_semantic_sim_onto`: Calculate the Resnik or Lin semantic similarity for an ontology.
- `calculate_wang_sim`: Calculate the Wang similarities between terms
- `calculate_wang_sim_mat`: Calculate the Wang similarity matrix
- `calculate_wang_sim_onto`: Calculate the Wang similarity for an ontology.
- `filter_similarities`: Filter the calculated similarities
- `calculate_critical_value`: Calculates the critical value
- `get_sim_matrix`: Get the similarity matrix
- `get_ontology_ancestry`: Return ancestry terms from an ontology

## Helpers for graph-based analysis

Different methods working on graphs; diffuse information over a network and identify communities, generate reciprocal best hit graphs from correlations or set similarities, fuse networks together via similarity network fusion.

- `NetworkDiffusions`: Network diffusion class
- `calculate_diffusion_auc`: Calculate the AUROC for a diffusion score
- `community_detection`: Identify privileged communities based on a given diffusion vector
- `constrained_page_rank`: Constrained personalised page rank
- `constrained_page_rank_ls`: Constrained personalised page rank over a list
- `diffuse_seed_nodes`: Diffuse seed genes over a network
- `permute_seed_nodes`: Generate permuation scores for the diffusion
- `get_diffusion_perms`: Get the diffusion permutations
- `generate_personalisation_vec`: Helper function to create personalisation vectors
- `tied_diffusion`: Diffuse seed genes in a tied manner over a network
- `RbhGraph`: Reciprocal best hit graph
- `find_rbh_communities`: Find RBH communities
- `generate_rbh_graph`: Generate an RBH graph.
- `get_diffusion_vector`: Get the diffusion vector
- `get_rbh_res`: Get the RBH results
- `SimilarityNetworkFusion`: Similarity network fusion
- `add_snf_data_modality`: Add a data modality for SNF generation
- `get_snf_params`: Get the SNF params
- `get_snf_final_mat`: Get the final SNF matrix
- `get_snf_adjcacency_mat`: Get an individual affinity matrix
- `params_graph_resolution`: Wrapper function to generate resolution parameters for Leiden or Louvain clustering.
- `run_snf`: Run the SNF algorithm
- `params_community_detection`: Wrapper function to generate community detection parameters
- `params_snf`: Wrapper function to generate SNF parameters

## Single cell class and getters

THE single cell class with a large number of getters.

- `SingleCells`: bixverse SingleCells class
- `SingleCellCountData`: Single cell count data handler
- `add_sc_new_obs`: Add an obs table derived from a method to the SingleCells.
- `get_sc_obs`: Getter the obs table
- `get_sc_var`: Getter the var table
- `get_sc_counts`: Getter the counts
- `get_available_embeddings`: Get the available embeddings
- `get_cell_indices`: Get the index position for a gene
- `get_cell_names`: Get the cell names
- `get_cells_to_keep`: Get the cells to keep
- `reset_cells_to_keep`: Reset the cells to keep
- `get_sc_cache_status`: Status of everything held in a single cell object's caches
- `check_sc_state`: Check that cached artefacts still match the object's state
- `assert_sc_state`: Assert that cached artefacts still match the object's state
- `get_embedding`: Get the embedding
- `get_gene_indices`: Get the index position for a gene
- `get_gene_names`: Get the gene names
- `get_cell_info`: Get the cell idx (R-based) and cell names
- `get_hvg`: Get the HVG
- `get_knn_mat`: Get the KNN matrix
- `get_knn_obj`: Get the KNN object
- `get_magic`: Get the MAGIC imputed layer
- `set_magic`: Set/add the MAGIC imputed layer
- `remove_magic`: Remove the MAGIC imputed layer
- `get_pca_singular_val`: Get the PCA singular values
- `get_pca_loadings`: Get the PCA loadings
- `get_pca_factors`: Get the PCA factors
- `get_snn_graph`: Get the sNN graph
- `get_gene_names_from_idx`: Get the gene names based on the gene idx
- `get_sc_available_features`: Returns the available features for single cell applications
- `setnames_sc`: Rename columns in the obs or var table
- `set_sc_new_obs_col`: Add a new column to the obs table
- `set_sc_new_obs_col_multiple`: Add multiple new columns to the obs table
- `set_sc_new_var_cols`: Add a new column to the var table
- `drop_cols_sc`: Drop columns from the obs or var table

## Single cell subsets and pipelines

Sub-clustering a group of cells and running the same chain per group.

- `SingleCellsSubset`: bixverse single cell subset class
- `merge_subset_obs`: Merge obs columns from subsets back into the parent object
- `sc_pipeline`: Construct an empty single cell pipeline
- `%>>%`: Append a step to a pipeline
- `apply_pipeline`: Apply a pipeline to a single cell object
- `apply_pipeline_per_group`: Apply a pipeline independently to each group of a SingleCells object
- `validate_pipeline`: Check that a pipeline can run on a given class
- `meta_cells_per_group`: Generate source-pure meta cells and merge them
- `step_hvg_sc`: Pipeline step: identify highly variable genes
- `step_pca_sc`: Pipeline step: PCA
- `step_neighbours_sc`: Pipeline step: nearest neighbours
- `step_clusters_sc`: Pipeline step: graph-based clustering
- `step_bbknn_sc`: Pipeline step: BBKNN batch correction
- `step_fast_mnn_sc`: Pipeline step: fastMNN batch correction
- `step_harmony_sc`: Pipeline step: Harmony batch correction
- `step_harmony_v2_sc`: Pipeline step: Harmony v2 batch correction
- `step_metacells_sc`: Pipeline step: generate meta cells

## Single cell class for multi-modal data and getters

This version allows you to also work ADT counts

- `SingleCellsMultiModal`: bixverse SingleCells (multi modal) class
- `new_adt_counts_clr`: Generates a new ADTCounts class
- `new_adt_counts_dsb`: Generates a new ADTCounts class via DSB normalisation
- `add_adt_counts_sc`: Add ADT counts to SingleCellsMultiModal
- `detect_adt_isotypes`: Detect likely isotype-control features by name pattern
- `get_adt_feature_info`: Get the ADT feature info
- `get_adt_names`: Get the ADT feature names
- `get_adt_sample_info`: Get the ADT sample info
- `params_sc_dsb`: Default parameters for DSB ADT normalisation
- `read_multi_tenx_h5_adt`: Read in 10x h5 ADT data from multiple files
- `read_tenx_h5_adt`: Read in 10x h5 ADT data
- `remove_adt_isotypes`: Return the ADT feature names removing the isotypes

## Meta cell-related generation, classes and methods

Generating meta cells.

- `MetaCells`: bixverse meta cell class
- `calc_diffusion_coordinates`: Calculate diffusion coordinates
- `calc_manifold_metrics`: Calculate manifold metrics
- `calc_meta_cell_purity`: Calculate meta cell purity
- `get_meta_cell_purity`: Calculate meta cell purity without mutating object state
- `generate_supercells_sc`: Generate SuperCells and return a MetaCells object
- `generate_bt_meta_cells_sc`: Generate meta cells based on hdWGCNA and return a MetaCells object
- `generate_seacells_sc`: Generate meta cells based on SEACells and return a MetaCells object
- `merge_meta_cells`: Merge meta cell objects into one
- `params_sc_supercell`: Wrapper function for parameters for SuperCell generation
- `params_sc_bt_metacells`: Wrapper function for parameters for bootstrapped meta cell generation
- `params_sc_seacells`: Wrapper function for the SEACells parameters

## Single cell i/o

I/O functions for single cell. Load in h5ad, mtx, Seurat or R data into Rust and the DuckDB supporting the metadata.

- `get_cell_ranger_params`: Helper to generate cell ranger input parameters
- `get_h5ad_dimensions`: Helper function to get the dimensions and storage format
- `prescan_h5ad_files`: Pre-scan multiple h5ad files for multi-sample loading
- `prescan_mtx_dirs`: Prescan multiple mtx directories for a multi-load
- `prescan_tenx_h5_files`: Pre-scan multiple 10x CellRanger h5 files for multi-sample loading
- `read_tenx_h5_metadata`: Read barcode and feature tables and metadata from a 10x h5 file
- `load_existing`: Load an existing SingleCells from disk
- `load_h5ad`: Load in h5ad to SingleCells
- `load_h5ad_norm`: Load in h5ad with normalised counts to SingleCells
- `load_mtx`: Load in mtx/plain text files to SingleCells
- `load_multi_mtx`: Load multiple mtx directories into a single SingleCells
- `load_multi_h5ad`: Load multiple h5ad files into a single SingleCells
- `stream_h5ad`: Stream in h5ad to SingleCells (alias)
- `load_r_data`: Load in data directly from R objects.
- `load_seurat`: Load in Seurat to SingleCells
- `load_tenx_h5`: Load in a 10x CellRanger h5 file to SingleCells
- `load_multi_tenx_h5`: Load multiple 10x CellRanger h5 files into a single SingleCells
- `read_h5ad_metadata`: Read obs and var tables and metadata from an h5ad file
- `read_h5ad_x_summary`: Read summary statistics from the X slot of an h5ad file
- `save_sc_exp_to_disk`: Save memory-bound data to disk
- `merge_sc_experiments`: Merge multiple SingleCells experiments into one
- `params_sc_min_quality`: Wrapper function to generate QC metric params for single cell
- `params_sc_mtx_io`: Wrapper function to provide data for mtx-based loading

## Single cell processing

Helpers to process single cell data. Doublet detection, proportions of gene sets, HVG (batch-aware), PCA and batch corrections.

- `scrublet_sc`: Doublet detection with Scrublet
- `call_doublets_manual`: Manually readjust Scrublet doublet call thresholds
- `doublet_detection_boost_sc`: Doublet detection with boosted doublet classification
- `scdblfinder_sc`: Run scDblFinder doublet detection on a SingleCells object
- `gene_set_proportions_sc`: Calculate the proportions of reads for specific gene sets
- `per_cell_qc_outlier`: Use MAD outlier detection on per-cell QC metrics
- `run_cell_qc`: Run outlier detection on per-cell QC metrics
- `run_cell_qc_fixed`: Fixed-threshold cell QC
- `rescue_cells`: Rescue MAD-flagged cells that fall within safe bounds
- `flag_cells`: Add hard-threshold flags to a CellQc object
- `find_hvg_sc`: Identify HVGs
- `find_hvg_batch_aware_sc`: Identify HVGs (batch aware)
- `get_hvg_data_sc`: Identify HVGs without mutating object state
- `calculate_pca_sc`: Run PCA for single cell
- `generate_sc_knn`: Generate a new SingleCellNearestNeighbour from data
- `find_neighbours_sc`: Find the neighbours for single cell.
- `run_magic_sc`: Impute a subset of genes with MAGIC
- `top_genes_perc_sc`: Calculate the proportions of reads for the Top N genes
- `params_sc_magic`: Wrapper function for MAGIC imputation parameters
- `params_norm_doublets_defaults`: Helper function to generate normalisation defaults for doublet detection.
- `params_boost`: Wrapper function for Boost parameters
- `params_sc_hvg`: Wrapper function for HVG detection parameters.
- `params_sc_pca`: Wrapper for PCA specifically designed for single cells
- `params_scrublet`: Wrapper function for Scrublet doublet detection parameters
- `params_sc_fast_cluster`: Fast single cell clustering parameters
- `params_sc_neighbours`: Wrapper function for parameters for neighbour identification in single cell
- `params_scdblfinder`: Wrapper function for scDblFinder doublet detection parameters
- `params_hvg_defaults`: Helper function to generate HVG defaults
- `params_pca_defaults`: Helper function to generate default parameters for PCA
- `params_knn_defaults`: Helper function to generate kNN defaults
- `params_sc_knn`: Parameters for single cell kNN searches
- `params_kmeans_defaults`: K-mean parameter defaults.
- `params_fast_cluster_default`: Helper function to generate default parameters for the fast clustering for the doublet detection methods

## Single cell batch correction methods

Batch correction methods and metrics for single cell

- `fast_mnn_sc`: Run fastMNN
- `harmony_sc`: Run Harmony
- `harmony_v2_sc`: Run Harmony v2
- `bbknn_sc`: Run BBKNN
- `seurat_cca_sc`: Run Seurat CCA integration
- `seurat_rpca_sc`: Run Seurat rPCA integration
- `calculate_kbet_sc`: Calculate kBET scores
- `calculate_batch_asw_sc`: Calculate batch average silhouette width
- `calculate_batch_lisi_sc`: Calculate batch LISI scores
- `params_sc_fastmnn`: Wrapper function for the fastMNN parameters
- `params_sc_harmony`: Default parameters for Harmony batch correction
- `params_sc_harmony_v2`: Default parameters for Harmony v2 batch correction
- `params_sc_bbknn`: Wrapper function for the BBKNN parameters
- `params_sc_seurat_cca`: Wrapper function for the Seurat CCA parameters
- `params_sc_seurat_rpca`: Wrapper function for the Seurat rPCA parameters

## Single cell analysis methods

A large number of different methods to extract insights from your single cell experiment. Gene set scoring, DGEs, kNN generations, pseudo-bulk count extraction, miloR, Hotspot, VISION and SCENIC.

- `aucell_sc`: Calculate AUC scores (akin to AUCell)
- `module_scores_sc`: Calculate module activity scores
- `fast_cluster_sc`: Run fast Louvain clustering on a SingleCells object
- `find_clusters_sc`: Graph-based clustering of cells on the sNN graph
- `find_markers_sc`: Calculate DGE between two cell groups
- `find_all_markers_sc`: Find all markers
- `find_specific_markers_sc`: Find markers that are specific to a cell group
- `get_pseudobulked_sc`: Generate pseudo-bulked matrices
- `generate_knn_sc`: Generate a SingleCellNearestNeighbour from a single cell class
- `get_differential_abundance_res`: Get the differential abundance results
- `hotspot_autocor_sc`: Calculate the local auto-correlation of a gene
- `hotspot_gene_cor_sc`: Calculate the local pairwise gene-gene correlation
- `generate_hotspot_membership`: Identify hotspot gene clusters
- `get_hotspot_membership`: Get the hotspot gene membership table
- `get_miloR_abundances_sc`: Generate an miloR abundance object for differential abundance testing
- `meld_sc`: Run MELD signal smoothing for differential abundance estimation
- `run_palantir_sc`: Run Palantir trajectory inference
- `run_paga_sc`: Run PAGA graph abstraction
- `run_gene_trends_sc`: Fit gene trends over Palantir pseudotime
- `get_index_cells`: Get the index cells
- `add_nhoods_info`: Add neighbourhood info on majority cell type
- `test_nhoods`: Test neighbourhoods for differential abundance
- `vision_sc`: Calculate VISION scores
- `vision_w_autocor_sc`: Calculate VISION scores (with auto-correlation scores)
- `identify_tf_to_genes`: Identify the TF to gene regulation
- `scenic_gene_filter_sc`: Filter genes for SCENIC GRN inference
- `scenic_grn_sc`: Run SCENIC GRN inference
- `get_cistarget_res`: Extract the TF to gene data from the ScenicGrn object
- `get_tf_to_gene`: Extract the TF to gene data from the ScenicGrn object
- `tf_to_genes_correlations`: Generate TF to gene correlations
- `tf_to_genes_motif_enrichment`: Run the SCENIC motif enrichment
- `binarise_regulon_activity`: Binarise regulon activity into on/off calls
- `build_regulons`: Build the final regulons
- `nmf_sc`: Run single-run NMF on single cell or meta cell data
- `stabilised_nmf_sc`: Run stabilised (multi-run) NMF on single cell or meta cell data
- `get_best_run`: Get the best run from a stabilised NMF result
- `get_w`: Get the W (gene loadings) matrix
- `get_h`: Get the H (cell activations) matrix
- `params_sc_aucell`: Wrapper function for parameters for AUCell
- `params_sc_hotspot`: Wrapper function for parameters for HotSpot
- `params_sc_miloR`: Wrapper function for parameters for MiloR
- `params_sc_vision`: Wrapper function for parameters for VISION with auto-correlation
- `params_scenic`: Constructor for SCENIC parameters
- `params_scenic_extra_trees_defaults`: Default parameters for the SCENIC ExtraTrees regression learner
- `params_scenic_gradient_boosting_defaults`: Default parameters for the SCENIC GradientBoosting (GRNBoost2) regression learner
- `params_scenic_random_forest_defaults`: Default parameters for the SCENIC RandomForest regression learner
- `params_scenic_binarise`: Wrapper function for parameters for the SCENIC binarisation
- `params_meld`: Constructor for MELD parameters
- `params_sc_palantir`: Wrapper function for Palantir parameters
- `params_sc_branch_selection`: Wrapper function for the branch cell selection parameters
- `params_sc_gene_trends`: Wrapper function for gene trend parameters
- `params_nmf_hals`: Wrapper function for NMF (HALS) parameters

## Single cell multi-modal analysis methods

Methods to analyse multi-modal single cell data

- `calculate_pca_adt_sc`: Calculate the PCA on top of the normalised ADT counts
- `generate_wnn_graph_sc`: Generate a weighted nearest neighbour (WNN) graph
- `params_sc_wnn`: Wrapper function for WNN parameters

## Single-cell related classes and methods

Additional helpers for specific small sub classes used in single cell.

- `calc_knn_metrics`: Calculate recall at k and distance ratio
- `get_centroids_sc`: Get k-means centroids from a fast cluster result
- `get_feature_mat`: Get the feature matrix used for the classifier
- `get_kmeans_clusters`: Get k-means cluster assignments from a fast cluster result
- `get_knn_dist`: Get the KNN distance
- `get_marker_summary`: Get the per-gene marker summaries across all rivals
- `get_marker_comparisons`: Get the per-rival marker statistics
- `get_data`: Get the ready obs data from various sub method
- `get_scores`: Get scores
- `new_sc_knn`: Helper function to generate kNN data with distances

## Reference mapping and cell type annotations in single cell

Helpers to do reference mapping and cell type identification in single cell

- `SymphonyReference`: bixverse SymphonyReference class
- `add_symphony_labels`: Add labels to a Symphony reference post-hoc
- `build_symphony_ref`: Build a Symphony reference from a SingleCells object
- `transfer_labels_symphony`: Transfer labels from a Symphony reference to a query via kNN majority vote
- `get_symphony_hvg_names`: Getter for the HVG gene names of a Symphony reference
- `get_symphony_labels`: Getter for the stored labels of a Symphony reference
- `get_symphony_loadings`: Getter for the PCA loadings of a Symphony reference
- `get_symphony_z_corr`: Getter for the corrected embedding of a Symphony reference
- `map_symphony_query`: Map a SingleCells query onto a Symphony reference
- `prepare_cell_markers`: Helper function to prepare cell markers
- `calc_sc_type_scores`: Calculate ScType scores per cell
- `score_clusters`: Score clusters based on ScType
- `assign_sc_type`: Assign cell types per cell based on ScType
- `params_sctype_cells`: Parameters for the per-cell ScType assignment
- `params_symphony_map`: Default parameters for Symphony query mapping

## Single cell ligand receptor

Classes, functions and generics/methods for ligand receptor analysis for single cell.

- `compute_expression_info_sc`: Compute per-cluster mean expression and expressing fraction for a gene set
- `generate_ligand_target_influence`: Generate the ligand to target influence matrix
- `get_influence`: Get the ligand-target influence matrix
- `ligand_activity_scores`: Compute ligand activity scores against gene sets
- `params_ligand_target`: Parameters for ligand to target influence computation
- `prioritise_interactions`: Prioritise sender-ligand-receiver-receptor interactions

## Single cell plotting stuff

Various helpers that generate data for plotting single cell, such as 2D embeddings, or extract specific columns from the metadata or genes (and their summaries) from the binary storage files.

- `sc_knn_to_nearest_neighbours`: Convert SingleCellNearestNeighbour to manifoldsR NearestNeighbours
- `umap_sc`: Run UMAP on a SingleCells/MetaCells object
- `tsne_sc`: Run t-SNE on a SingleCells/MetaCells object
- `phate_sc`: Run PHATE on a SingleCells/MetaCells object
- `extract_dot_plot_data`: Extract grouped gene statistics for dot plots
- `extract_gene_expression`: Extract normalised gene expression for plotting
- `extract_embedding_data`: Extract embedding coordinates for plotting
- `extract_feature_pair`: Extract a pair of features for scatter / hex plots
- `extract_feature_plot_data`: Extract per-cell expression mapped onto an embedding
- `extract_gene_violin_data`: Extract per-cell expression grouped for violin plots
- `extract_paga_plot_data`: Extract the PAGA graph positioned on an embedding

## Statistical functions

Any types of functions that help with statistics

- `calculate_effect_size`: Calculate the Hedge G effect between two matrices
- `calculate_tom`: Calculate the TOM from a correlation matrix
- `calculate_tom_from_exp`: Calculate the TOM from an expression matrix
- `f1_score_confusion_mat`: F1 scores on top of a confusion matrix
- `fast_ica_rust`: Fast ICA via Rust
- `fast_ica_rust_helper`: Fast ICA via Rust from processed data
- `get_inflection_point`: Identify the inflection point for elbow-like data
- `ot_harmonic_score`: Calculates a harmonic sum normalised between 0 to 1.
- `robust_scale`: Robust scaler.

## Plotting helpers

Some core plotting helpers in the package (usually for QC). The ones to plot downstream results can be found in bixverse.plots.

- `plot_boxplot_normalization`: Helper plot function for boxplot of normalized data
- `plot_epsilon_res`: Plot the epsilon vs. power law goodness of fit result
- `plot_hvgs`: Plot the highly variable genes
- `plot_ica_ncomp_params`: Plot various parameters with no comp
- `plot_ica_stability_individual`: Plot the stability of the ICA components
- `plot_optimal_cuts`: Plot the k cuts vs median R2
- `plot_pca`: Helper plot function for pca with contrasts
- `plot_pca_res`: Plot the PCA data
- `plot_preprocessing_genes`: Helper plot function of distribution of genes by samples
- `plot_preprocessing_outliers`: Helper plot function for identification of outliers
- `plot_rbf_impact`: Helper function to plot distance to affinity relationship
- `plot_resolution_res`: Plot the resolution results.
- `plot_voom_normalization`: Helper plot function for Voom normalisation

## Data downloads and synthetic data generation

Functions and helpers to download or generate synthetic data.

- `download_cd34_data`: Download the CD34 example data from SEACells
- `download_marrow_cd34`: Download the marrow CD34 example data from Palantir
- `download_pbmc3k`: Download PBMC3K data from Zenodo
- `download_demuxlet_pbmc`: Download PBMCs with demuxlet doublet information
- `download_pbmc_batches`: Download two different PBMC data sets for batch correction testing
- `download_pbmc_totalseq_data`: Download the PBMC TotalSeq data with ADT counts
- `download_pbmc8k`: Download PBMC8K data from Zenodo
- `calculate_sparsity_stats`: Helper function to calculate the induced sparsity
- `generate_gene_module_data`: Generates synthetic gene module data.
- `generate_single_cell_test_data`: Single cell test data
- `cell_cycle_genes`: Cell cycle genes
- `write_cellranger_output`: Helper function to write data to a cell ranger like output
- `write_h5ad_sc`: Helper function to write data to h5ad format
- `write_h5ad_sc_dense`: Helper function to write data to a dense h5ad file
- `params_sc_synthetic_data`: Default parameters for generation of synthetic single cell data (RNA)
- `params_synthetic_bulk_rnaseq`: Wrapper function to generate synthetic bulk RNAseq parameters
- `params_bulk_sparsity`: Wrapper function to generate bulk sparsification parameters
- `synthetic_signal_matrix`: Generates a simple synthetic, pseudo gene expression matrix
- `simulate_dropouts`: Simulate sequencing-depth dropouts on synthetic bulk data
- `synthetic_bulk_cor_matrix`: Generates synthetic bulk RNAseq data
- `synthetic_c_pca_data`: Generates synthetic data for contrastive PCA exploration.

## Utils

All types of other random helpers without a clear pattern

- `AnnDataParser`: Class for Anndata
- `calculate_sparsity_stats`: Helper function to calculate the induced sparsity
- `find_threshold_otsu`: Find a threshold via the Otsu method
- `get_seurat_counts_to_list`: Transform Seurat raw counts into a List
- `install_agent_skill`: Install the bixverse agent skill
- `knn_graph_label_propagation`: kNN-based graph label propagation
- `params_label_propagation`: Wrapper function to generate label propagation parameters
- `upper_triangle_to_sparse`: Transform an upper triangle-stored matrix to a sparse one
- `upper_triangular_sym_mat`: Class for symmetric correlation matrices
- `to_snake_case`: Helper function to transform strings to snake_case

## Rust wrappers

Everything Rusty - only use this if you know what you are doing... Maybe useful for your own package? Use with care and read the documentation! The ones exposed here are general enough to be useful in other packages. There is a lot more under the hood...

92 `rs_*` functions are exposed here. They are the raw extendr bindings with no input validation. Use the R wrapper instead; only reach for these if you are building on top of bixverse and know exactly what you are doing.

## Not on the package website

Exported but absent from `_pkgdown.yml`. Mostly internal constructors and `rs_*` bindings that take on-disk streaming input. Usable, but undocumented on the website, so read the roxygen with `?fn` first.

- `BixverseBaseClass`: BixverseBaseClass
- `bulk_coexp`: Bulk RNAseq co-expression modules (deprecated)
- `bulk_dge`: Bulk RNAseq differential gene expression class (deprecated)
- `detect_raw_count_slot`: Detect which slot holds raw integer counts in an h5ad file
- `gene_ontology_data`: Gene Ontology data (deprecated)
- `generate_single_cell_test_data_adt`: Single cell test data (ADT)
- `get_params.Hotspot`: Get the parameters that were used.
- `get_params.miloR`: Get the parameters that were used.
- `get_params.NmfResult`: Get the parameters that were used.
- `get_params.ScenicGrn`: Get the parameters that were used.
- `get_params.StabilisedNmfResult`: Get the parameters that were used.
- `get_sc_cache`: Getter the memory-stored data from the class
- `get_sc_duckdb`: Getter for the single cell DuckDB connection
- `get_sc_map`: Getter for the different maps in the object
- `get_sc_rust_ptr`: Getter for the single cell Rust pointer
- `mc_counts_to_list`: Transform the counts to a Rust-specific list
- `mc_get_clr_offsets`: Get the offsets for the CLR/PFlogPF transformation prior PCA
- `meta_cells`: bixverse meta cell class (deprecated)
- `network_diffusions`: Network diffusion class (deprecated)
- `new_gene_trends_res`: Helper function to generate the gene trend results
- `new_nmf_result`: Constructor for single-run NMF results
- `new_paga_res`: Helper function to generate the PAGA results
- `new_palantir_res`: Helper function to generate the Palantir results
- `new_sc_cache`: Helper function to hold relevant cached data
- `new_sc_hotspot_res`: Helper function to generate HotSpot data
- `new_sc_list`: Generate an ScListRes instance
- `new_sc_magic`: Helper function to generate the MAGIC imputed layer
- `new_sc_mapper`: Helper function to construct relevant maps
- `new_sc_matrix`: Generate an ScMatrixRes instance
- `new_sc_specific_markers`: Constructor for the one-vs-many specific marker results
- `new_scenic_grn`: Constructor for SCENIC GRN results
- `new_stabilised_nmf_result`: Constructor for stabilised (multi-run) NMF results
- `ontology`: Ontology class (deprecated)
- `params_sc_synthetic_data_adt`: Default parameters for generation of synthetic single cell data (ADT)
- `prep_data_gower_hamming_dist`: Transform data.tables into matrices for distance calculations
- `rbh_graph`: Reciprocal best hit graph (deprecated)
- `remove_knn`: Remove the KNN data
- `remove_snn_graph`: Remove the sNN graph
- `rs_adt_clr`: Applies CLR normalisation on ADT counts
- `rs_aucell`: Calculate AUCell in Rust
- `rs_bbknn`: BBKNN implementation in Rust
- `rs_bbknn_filtering`: Reduce BBKNN results to Top X neighbours
- `rs_build_symphony_ref`: Build a Symphony reference (Rust)
- `rs_calc_es`: Calculates the traditional GSEA enrichment score
- `rs_calc_gsea_stat_cumulative_batch`: Helper function to generate fgsea simple-based permutations
- `rs_calc_gsea_stat_traditional_batch`: Helper function to generate traditional GSEA-based permutations
- `rs_calc_gsea_stats`: Rust implementation of the fgsea::calcGseaStat() function
- `rs_calc_multi_level`: Calculates p-values for pre-processed data
- `rs_calculate_dge_mann_whitney`: Calculate DGEs between cells based on Mann Whitney stats
- `rs_calculate_dge_one_vs_many`: Calculate one-vs-many AUROC DGEs for specific markers
- `rs_compute_cluster_expr_stats`: Compute cluster statistics for NicheNet prioritisation
- `rs_create_random_aucs`: Create random AUCs
- `rs_critval`: Calculate the critical value
- `rs_critval_mat`: Calculate the critical value
- `rs_dsb`: Run DSB normalisation on raw ADT counts
- `rs_extract_counts_plots`: Helper to extract single cell counts as a dense vector for plotting
- `rs_extract_grouped_gene_stats`: Calculates the gene statistics for a set of cell groups and genes
- `rs_extract_several_genes_plots`: Helper to extract single cell counts for several genes
- `rs_fast_cluster_sc`: Runs fast Louvain cluster on the data
- `rs_fast_cluster_sc_grid`: Runs fast Louvain cluster on the data (with multiple seeds)
- `rs_filter_onto_sim`: Filter the term similarities for a specific critical value
- `rs_generate_bulk_rnaseq`: Generation of bulkRNAseq-like data with optional correlation structure
- `rs_get_gs_indices`: Helper function to rapidly retrieve the indices of the gene set members
- `rs_get_metacells_bootstrapped`: Generate meta cells (hdWGCNA method)
- `rs_get_seacells`: Generate SEACells
- `rs_hotspot_autocor`: Calculate gene spatial auto-correlations
- `rs_hotspot_cluster_genes`: Cluster the genes by Z-score together
- `rs_hotspot_gene_cor`: Calculate gene to gene spatial correlations
- `rs_importance_threshold`: SCENIC: Select TF-gene pairs by per-gene importance threshold
- `rs_knn_mat_to_edge_pairs`: Flatten kNN matrix to edge list
- `rs_magic_impute`: Impute a subset of genes with MAGIC
- `rs_make_milor_nhoods`: Generate the neighbourhoods akin to the miloR approach
- `rs_meld_sc`: Run MELD
- `rs_metacell_compactness`: Calculates the compactness of the MetaCells based on diffusion map coordinates
- `rs_metacell_density`: Calculates diffusion maps for density calculations for meta cells
- `rs_metacell_separation`: Calculates the separation of the centroids of the MetaCells based on diffusion map coordinates.
- `rs_mnn`: FastMNN batch correction in Rust
- `rs_module_scoring`: Calculate module activity scores in Rust
- `rs_nmf_multi_bulk`: Run multiple NMF (HALS) restarts on a bulk expression matrix
- `rs_nmf_multi_sc`: Run multiple NMF (HALS) restarts over a set of single cells and genes
- `rs_nmf_single_bulk`: Run NMF (HALS) on a bulk expression matrix (single run)
- `rs_nmf_single_sc`: Run NMF (HALS) over a set of single cells and genes
- `rs_paga`: Run PAGA over a kNN graph and a clustering
- `rs_pairwise_gene_cors`: Calculates pairwise gene correlations in single cell
- `rs_pairwise_gene_cors_mc`: Calculate the pairwise gene-correlation for meta cells
- `rs_palantir`: Run Palantir trajectory inference over a kNN graph
- `rs_prepare_gsva_gs`: Prepare a pathway list for GSVA
- `rs_pseudobulk_cells_dense`: Pseudo-bulk a set of cells (dense)
- `rs_pseudobulk_cells_sparse`: Pseudo-bulk a set of cells (sparse)
- `rs_rbf_iterate_epsilons`: Helper to identify the right epsilon parameter
- `rs_read_tenx_h5_modality`: Loads in a modality from a 10x h5 file
- `rs_regulon_thresholds`: Derive the on/off threshold per regulon in Rust
- `rs_sample_ids_for_cell_types`: Helper function to generate sample identifiers based on cells
- `rs_sc_doublet_detection`: Detect Doublets via BoostClassifier (in Rust)
- `rs_sc_get_gene_set_perc`: Calculate the percentage of gene sets in the cells
- `rs_sc_get_top_genes_perc`: Calculates the cumulative proportion of the top X genes
- `rs_sc_hvg`: Calculate the percentage of gene sets in the cells
- `rs_sc_hvg_batch_aware`: Calculate HVG per batch
- `rs_sc_otsu_method`: Run Otsu's method
- `rs_sc_pca`: Calculates PCA for single cell
- `rs_sc_pca_sparse`: Calculates sparse PCA for single cell
- `rs_sc_scdblfinder`: Run scDblFinder doublet detection
- `rs_sc_scrublet`: Scrublet Rust interface
- `rs_sc_type`: Run the ScType scoring approach
- `rs_sc_type_assign_cells`: Per-cell ScType assignment with optional kNN smoothing
- `rs_sc_type_cluster_assignment`: Score the individual clusters based on ScType
- `rs_scenic_gene_filter`: Identifies genes to include into a SCENIC analysis
- `rs_scenic_grn`: SCENIC: Generating gene-regulatory networks
- `rs_scenic_grn_streaming`: SCENIC: Generating gene-regulatory networks (streaming version)
- `rs_seurat_cca`: Seurat CCA batch correction in Rust
- `rs_seurat_rpca`: Seurat rPCA batch correction in Rust
- `rs_simple_and_multi_err`: Calculates the simple and multi error for fgsea multi level
- `rs_simulate_dropouts`: Sparsify bulkRNAseq like data
- `rs_snf_affinity_cat`: Calculate the SNF affinity matrix for categorical values
- `rs_snf_affinity_continuous`: Calculate the SNF affinity matrix for continuous values
- `rs_snf_affinity_mixed`: Calculate the SNF affinity matrix for mixed values
- `rs_split_cor_signs`: Helper function to split correlation matrices by sign
- `rs_supercell`: Generate SuperCells.
- `rs_symphony_map_query`: Map a query onto a Symphony reference (Rust)
- `rs_synthetic_sc_adt_with_cell_types`: Generates synthetic ADT counts with defined cell types
- `rs_synthetic_sc_data_with_cell_types`: Generates synthetic data for single cell
- `rs_top_k_targets`: SCENIC: Select the Top TF <> Gene pairs
- `rs_transfer_labels_symphony`: Transfer labels from a Symphony reference to a query via kNN vote
- `rs_mc_vision`: Calculate VISION pathway scores in Rust (for meta cells)
- `rs_mc_vision_with_autocorrelation`: Calculate VISION pathway scores in Rust with auto-correlation (for meta cells)
- `rs_vision`: Calculate VISION pathway scores in Rust
- `rs_vision_with_autocorrelation`: Calculate VISION pathway scores in Rust with auto-correlation
- `set_cell_mapping`: Set cell mapping
- `set_cells_to_keep`: Set cells to keep
- `set_embedding`: Add additional embeddings to the class
- `set_gene_mapping`: Set gene mapping
- `set_hvg`: Set the HVG genes
- `set_knn`: Set/add KNN
- `set_pca_factors`: Set/add PCA factors
- `set_pca_loadings`: Set/add PCA loadings
- `set_pca_singular_vals`: Set/add PCA singular values
- `set_snn_graph`: Set/add KNN
- `single_cell_exp`: bixverse single cell class (deprecated)
- `SingleCellDuckDB`: Class for storing single cell experimental data in DuckDB (nightly!)
- `SingleCellDuckDBBase`: Base class for the single cell DuckDB connection
- `snf`: Similarity network fusion (deprecated)
- `sparse_list_to_mat`: Helper function to transform the Rust-exported sparse matrices into R ones
- `testSNFParams`: Check SNF parameters
- `write_tenx_h5_sc`: Helper function to write data to a 10x CellRanger-style h5 file

