# Package index

## General getters

All types of general getters that work across various classes related to
co-expression module detection, graph-based clustering, etc.

- [`get_metadata()`](https://gregorlueg.github.io/bixverse/reference/get_metadata.md)
  : Return the metadata
- [`get_params()`](https://gregorlueg.github.io/bixverse/reference/get_params.md)
  [`get_params.Hotspot()`](https://gregorlueg.github.io/bixverse/reference/get_params.md)
  [`get_params.miloR()`](https://gregorlueg.github.io/bixverse/reference/get_params.md)
  [`get_params.ScenicGrn()`](https://gregorlueg.github.io/bixverse/reference/get_params.md)
  [`get_params.NmfResult()`](https://gregorlueg.github.io/bixverse/reference/get_params.md)
  [`get_params.StabilisedNmfResult()`](https://gregorlueg.github.io/bixverse/reference/get_params.md)
  : Get the parameters that were used.
- [`get_results()`](https://gregorlueg.github.io/bixverse/reference/get_results.md)
  : Get the final results from the class

## Gene set enrichment helpers

Everything and anything you need to do various types of gene set
enrichments; hypergeometric tests, GSVA, (ss)GSEA.

- [`gse_hypergeometric()`](https://gregorlueg.github.io/bixverse/reference/gse_hypergeometric.md)
  : Gene set enrichment (GSE) based on a hypergeometric test.
- [`gse_hypergeometric_list()`](https://gregorlueg.github.io/bixverse/reference/gse_hypergeometric_list.md)
  : Gene set enrichment (GSE) based on a hypergeometric test over a
  list.
- [`calc_fgsea()`](https://gregorlueg.github.io/bixverse/reference/calc_fgsea.md)
  : Bixverse implementation of the fgsea algorithm
- [`calc_fgsea_simple()`](https://gregorlueg.github.io/bixverse/reference/calc_fgsea_simple.md)
  : Bixverse implementation of the simple fgsea algorithm
- [`calc_gsea_traditional()`](https://gregorlueg.github.io/bixverse/reference/calc_gsea_traditional.md)
  : Bixverse implementation of the traditional GSEA algorithm
- [`calc_mitch()`](https://gregorlueg.github.io/bixverse/reference/calc_mitch.md)
  : Calculate a mitch gene set enrichments on contrast
- [`calc_gsva()`](https://gregorlueg.github.io/bixverse/reference/calc_gsva.md)
  : Bixverse implementation of GSVA
- [`calc_ssgsea()`](https://gregorlueg.github.io/bixverse/reference/calc_ssgsea.md)
  : Bixverse implementation of ssGSEA
- [`calc_singscore()`](https://gregorlueg.github.io/bixverse/reference/calc_singscore.md)
  : Bixverse implementation of singscore (single gene set)
- [`calc_singscore_multi()`](https://gregorlueg.github.io/bixverse/reference/calc_singscore_multi.md)
  : Bixverse implementation of singscore (multiple gene sets)
- [`calc_singscore_rank()`](https://gregorlueg.github.io/bixverse/reference/calc_singscore_rank.md)
  : Rank an expression matrix for singscore
- [`params_gsea()`](https://gregorlueg.github.io/bixverse/reference/params_gsea.md)
  : Wrapper function to generate GSEA parameters
- [`params_gsva()`](https://gregorlueg.github.io/bixverse/reference/params_gsva.md)
  : Wrapper function to generate GSVA parameters
- [`params_ssgsea()`](https://gregorlueg.github.io/bixverse/reference/params_ssgsea.md)
  : Wrapper function to generate ssGSEA parameters

## Gene Ontology-related gene set enrichment helpers

GSE methods specifically designed with the DAG of the Gene Ontology in
mind, identifying the most relevant Gene Ontology terms.

- [`load_go_human_data()`](https://gregorlueg.github.io/bixverse/reference/load_go_human_data.md)
  : Get the Gene Ontology data human
- [`get_go_data_human()`](https://gregorlueg.github.io/bixverse/reference/get_go_data_human.md)
  : Wrapper function to load and process the gene ontology data.
- [`process_go_data()`](https://gregorlueg.github.io/bixverse/reference/process_go_data.md)
  : Process Gene Ontology data into the right format
- [`GeneOntologyElim()`](https://gregorlueg.github.io/bixverse/reference/GeneOntologyElim.md)
  : Gene Ontology data
- [`gse_go_elim_method()`](https://gregorlueg.github.io/bixverse/reference/gse_go_elim_method.md)
  : Run gene ontology enrichment with elimination method.
- [`gse_go_elim_method_list()`](https://gregorlueg.github.io/bixverse/reference/gse_go_elim_method_list.md)
  : Run gene ontology enrichment with elimination method over a list.
- [`fgsea_go_elim()`](https://gregorlueg.github.io/bixverse/reference/fgsea_go_elim.md)
  : Run GO enrichment with elimination method over a continuous vectors
- [`fgsea_simple_go_elim()`](https://gregorlueg.github.io/bixverse/reference/fgsea_simple_go_elim.md)
  : Run GO enrichment with elimination with fgsea simple
- [`simplify_hypergeom_res()`](https://gregorlueg.github.io/bixverse/reference/simplify_hypergeom_res.md)
  : Simplify gene set results via ontologies

## CisTarget-related helpers

Functions related to CisTarget enrichment.

- [`download_cistarget_hg38()`](https://gregorlueg.github.io/bixverse/reference/download_cistarget_hg38.md)
  : Download CisTarget reference files for human (hg38)
- [`read_motif_annotation_file()`](https://gregorlueg.github.io/bixverse/reference/read_motif_annotation_file.md)
  : Read in the motif annotation file
- [`read_motif_ranking()`](https://gregorlueg.github.io/bixverse/reference/read_motif_ranking.md)
  : Read in the motif rankings and transform them into a matrix
- [`run_cistarget()`](https://gregorlueg.github.io/bixverse/reference/run_cistarget.md)
  : Main function to run CisTarget
- [`params_cistarget()`](https://gregorlueg.github.io/bixverse/reference/params_cistarget.md)
  : Wrapper function to CisTarget parameters

## Correlation-based gene module detections

Methods to identify co-expression modules via correlations in the data.
Includes methods for differential correlations, some graph-based and
some hierarchical clustering based ones.

- [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md)
  : Bulk RNAseq co-expression modules
- [`preprocess_bulk_coexp()`](https://gregorlueg.github.io/bixverse/reference/preprocess_bulk_coexp.md)
  : Process the raw data
- [`cor_module_check_epsilon()`](https://gregorlueg.github.io/bixverse/reference/cor_module_check_epsilon.md)
  : Iterate through different epsilon parameters
- [`cor_module_coremo_clustering()`](https://gregorlueg.github.io/bixverse/reference/cor_module_coremo_clustering.md)
  : Generates CoReMo-based gene modules
- [`cor_module_coremo_cor_sign()`](https://gregorlueg.github.io/bixverse/reference/cor_module_coremo_cor_sign.md)
  : Split CoReMo modules by correlation sign
- [`cor_module_coremo_eigengene()`](https://gregorlueg.github.io/bixverse/reference/cor_module_coremo_eigengene.md)
  : Calculate Eigengenes for CoReMo modules
- [`cor_module_coremo_stability()`](https://gregorlueg.github.io/bixverse/reference/cor_module_coremo_stability.md)
  : Assesses CoReMo-based gene module stability
- [`cor_module_graph_check_res()`](https://gregorlueg.github.io/bixverse/reference/cor_module_graph_check_res.md)
  : Iterate through Leiden resolutions for graph-based community
  detection.
- [`cor_module_graph_final_modules()`](https://gregorlueg.github.io/bixverse/reference/cor_module_graph_final_modules.md)
  : Identify correlation-based gene modules via graphs
- [`cor_module_processing()`](https://gregorlueg.github.io/bixverse/reference/cor_module_processing.md)
  : Prepare correlation-based module detection
- [`cor_module_tom()`](https://gregorlueg.github.io/bixverse/reference/cor_module_tom.md)
  : Update the correlation matrix to a TOM
- [`diffcor_module_processing()`](https://gregorlueg.github.io/bixverse/reference/diffcor_module_processing.md)
  : Prepare differential correlation-based module detection
- [`get_cor_graph()`](https://gregorlueg.github.io/bixverse/reference/get_cor_graph.md)
  : Get correlation-based graph
- [`get_diffcor_graph()`](https://gregorlueg.github.io/bixverse/reference/get_diffcor_graph.md)
  : Get differential correlation-based graph
- [`get_epsilon_res()`](https://gregorlueg.github.io/bixverse/reference/get_epsilon_res.md)
  : Return the epsilon data
- [`get_resolution_res()`](https://gregorlueg.github.io/bixverse/reference/get_resolution_res.md)
  : Return the resolution results
- [`get_outputs()`](https://gregorlueg.github.io/bixverse/reference/get_outputs.md)
  : Return the outputs
- [`params_cor_graph()`](https://gregorlueg.github.io/bixverse/reference/params_cor_graph.md)
  : Wrapper function for graph generation
- [`params_coremo()`](https://gregorlueg.github.io/bixverse/reference/params_coremo.md)
  : Wrapper function to generate CoReMo parameters

## Matrix factorisation methods

Methods to identify co-expression via matrix factorisations. Uses the
same class as the correlation-based ones. You have contrastive PCA, ICA
and dual graph-regularised dictionary learning as methods.

- [`contrastive_pca_processing()`](https://gregorlueg.github.io/bixverse/reference/contrastive_pca_processing.md)
  : Prepare class for contrastive PCA
- [`c_pca_plot_alphas()`](https://gregorlueg.github.io/bixverse/reference/c_pca_plot_alphas.md)
  : Plot various alphas for the contrastive PCA
- [`contrastive_pca()`](https://gregorlueg.github.io/bixverse/reference/contrastive_pca.md)
  : Apply contrastive PCA.
- [`get_c_pca_factors()`](https://gregorlueg.github.io/bixverse/reference/get_c_pca_factors.md)
  : Get the contrastive PCA factors
- [`get_c_pca_loadings()`](https://gregorlueg.github.io/bixverse/reference/get_c_pca_loadings.md)
  : Get the contrastive PCA loadings
- [`dgrdl_grid_search()`](https://gregorlueg.github.io/bixverse/reference/dgrdl_grid_search.md)
  : Grid search over DGRDL parameters
- [`dgrdl_result()`](https://gregorlueg.github.io/bixverse/reference/dgrdl_result.md)
  : Run DGRDL with the specified parameters
- [`get_ica_stability_res()`](https://gregorlueg.github.io/bixverse/reference/get_ica_stability_res.md)
  : Get the ICA component data (stability, convergence, nMI)
- [`get_grid_search_res()`](https://gregorlueg.github.io/bixverse/reference/get_grid_search_res.md)
  : Get the grid search results
- [`ica_evaluate_comp()`](https://gregorlueg.github.io/bixverse/reference/ica_evaluate_comp.md)
  : Iterate over different ncomp parameters for ICA
- [`ica_optimal_ncomp()`](https://gregorlueg.github.io/bixverse/reference/ica_optimal_ncomp.md)
  : Identify stability inflection point
- [`ica_processing()`](https://gregorlueg.github.io/bixverse/reference/ica_processing.md)
  : Prepare class for ICA
- [`ica_stabilised_results()`](https://gregorlueg.github.io/bixverse/reference/ica_stabilised_results.md)
  : Run stabilised ICA with a given number of components
- [`params_dgrdl()`](https://gregorlueg.github.io/bixverse/reference/params_dgrdl.md)
  : Wrapper function to generate DGRDL parameters
- [`params_ica_general()`](https://gregorlueg.github.io/bixverse/reference/params_ica_general.md)
  : Wrapper function for standard ICA parameters
- [`params_ica_ncomp()`](https://gregorlueg.github.io/bixverse/reference/params_ica_ncomp.md)
  : Wrapper function for ICA ncomp iterations
- [`params_ica_randomisation()`](https://gregorlueg.github.io/bixverse/reference/params_ica_randomisation.md)
  : Wrapper function for ICA randomisation

## Helpers for differential gene expression

Methods to help out with differential gene expression analyses in a
structured way. Useful when you have to analyse 10’s to 100’s of
differential gene expression results.

- [`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md)
  : Bulk RNAseq differential gene expression class
- [`add_new_metadata()`](https://gregorlueg.github.io/bixverse/reference/add_new_metadata.md)
  : Replace the meta data
- [`change_gene_identifier()`](https://gregorlueg.github.io/bixverse/reference/change_gene_identifier.md)
  : Change the primary gene identifier of BulkDge
- [`update_metadata_values()`](https://gregorlueg.github.io/bixverse/reference/update_metadata_values.md)
  : Replace values in a metadata column
- [`fix_meta_data_column()`](https://gregorlueg.github.io/bixverse/reference/fix_meta_data_column.md)
  : Helper to fix meta-data columns to be R conform
- [`remove_samples()`](https://gregorlueg.github.io/bixverse/reference/remove_samples.md)
  : Remove samples from object
- [`qc_bulk_dge()`](https://gregorlueg.github.io/bixverse/reference/qc_bulk_dge.md)
  : QC on the bulk dge data
- [`preprocess_bulk_dge()`](https://gregorlueg.github.io/bixverse/reference/preprocess_bulk_dge.md)
  : QC on the bulk dge data (DEPRECATED!)
- [`normalise_bulk_dge()`](https://gregorlueg.github.io/bixverse/reference/normalise_bulk_dge.md)
  : Normalise the count data for DGE.
- [`batch_correction_bulk_dge()`](https://gregorlueg.github.io/bixverse/reference/batch_correction_bulk_dge.md)
  : Run a linear batch correction
- [`bulk_dge_from_h5ad()`](https://gregorlueg.github.io/bixverse/reference/bulk_dge_from_h5ad.md)
  : Wrapper function to generate BulkDge object from h5ad
- [`calculate_dge_hedges()`](https://gregorlueg.github.io/bixverse/reference/calculate_dge_hedges.md)
  : Calculates the Hedge's G effect size
- [`calculate_all_dges()`](https://gregorlueg.github.io/bixverse/reference/calculate_all_dges.md)
  : Calculate all possible DGE variants (DEPRECATED!)
- [`calculate_dge_limma()`](https://gregorlueg.github.io/bixverse/reference/calculate_dge_limma.md)
  : Calculates the Limma Voom DGE
- [`calculate_pca_bulk_dge()`](https://gregorlueg.github.io/bixverse/reference/calculate_pca_bulk_dge.md)
  : Calculate PCA on the expression.
- [`calculate_rpkm()`](https://gregorlueg.github.io/bixverse/reference/calculate_rpkm.md)
  : RPKM calculation
- [`calculate_tpm()`](https://gregorlueg.github.io/bixverse/reference/calculate_tpm.md)
  : TPM calculation
- [`run_limma_voom()`](https://gregorlueg.github.io/bixverse/reference/run_limma_voom.md)
  : Wrapper for a Limma Voom analysis
- [`hedges_g_dge()`](https://gregorlueg.github.io/bixverse/reference/hedges_g_dge.md)
  : Calculate the effect size
- [`get_dge_effect_sizes()`](https://gregorlueg.github.io/bixverse/reference/get_dge_effect_sizes.md)
  : Return the effect size results
- [`get_dge_limma_voom()`](https://gregorlueg.github.io/bixverse/reference/get_dge_limma_voom.md)
  : Return the Limma Voom results
- [`get_dge_list()`](https://gregorlueg.github.io/bixverse/reference/get_dge_list.md)
  : Return the DGEList
- [`get_dge_qc_plot()`](https://gregorlueg.github.io/bixverse/reference/get_dge_qc_plot.md)
  : Return QC plots
- [`get_fpkm_counts()`](https://gregorlueg.github.io/bixverse/reference/get_fpkm_counts.md)
  : Return the FPKM-normalised counts
- [`get_gene_lengths()`](https://gregorlueg.github.io/bixverse/reference/get_gene_lengths.md)
  : Get the gene lengths
- [`get_model_fit()`](https://gregorlueg.github.io/bixverse/reference/get_model_fit.md)
  : Get the fitted model
- [`get_tpm_counts()`](https://gregorlueg.github.io/bixverse/reference/get_tpm_counts.md)
  : Return the TPM-normalised counts

## Biomedical ontologies

For dealing with ontologies and calculating (semantic) similarities in
disease, phenotype or gene ontologies.

- [`OntologySim()`](https://gregorlueg.github.io/bixverse/reference/OntologySim.md)
  : OntologySim class
- [`pre_process_sim_onto()`](https://gregorlueg.github.io/bixverse/reference/pre_process_sim_onto.md)
  : Pre-process data for subsequent ontology similarity
- [`calculate_information_content()`](https://gregorlueg.github.io/bixverse/reference/calculate_information_content.md)
  : Calculate the information content for each ontology term
- [`calculate_semantic_sim()`](https://gregorlueg.github.io/bixverse/reference/calculate_semantic_sim.md)
  : Calculate the Resnik or Lin semantic similarity
- [`calculate_semantic_sim_mat()`](https://gregorlueg.github.io/bixverse/reference/calculate_semantic_sim_mat.md)
  : Calculate the Resnik or Lin semantic similarity matrix
- [`calculate_semantic_sim_onto()`](https://gregorlueg.github.io/bixverse/reference/calculate_semantic_sim_onto.md)
  : Calculate the Resnik or Lin semantic similarity for an ontology.
- [`calculate_wang_sim()`](https://gregorlueg.github.io/bixverse/reference/calculate_wang_sim.md)
  : Calculate the Wang similarities between terms
- [`calculate_wang_sim_mat()`](https://gregorlueg.github.io/bixverse/reference/calculate_wang_sim_mat.md)
  : Calculate the Wang similarity matrix
- [`calculate_wang_sim_onto()`](https://gregorlueg.github.io/bixverse/reference/calculate_wang_sim_onto.md)
  : Calculate the Wang similarity for an ontology.
- [`filter_similarities()`](https://gregorlueg.github.io/bixverse/reference/filter_similarities.md)
  : Filter the calculated similarities
- [`calculate_critical_value()`](https://gregorlueg.github.io/bixverse/reference/calculate_critical_value.md)
  : Calculates the critical value
- [`get_sim_matrix()`](https://gregorlueg.github.io/bixverse/reference/get_sim_matrix.md)
  : Get the similarity matrix
- [`get_ontology_ancestry()`](https://gregorlueg.github.io/bixverse/reference/get_ontology_ancestry.md)
  : Return ancestry terms from an ontology

## Helpers for graph-based analysis

Different methods working on graphs; diffuse information over a network
and identify communities, generate reciprocal best hit graphs from
correlations or set similarities, fuse networks together via similarity
network fusion.

- [`NetworkDiffusions()`](https://gregorlueg.github.io/bixverse/reference/NetworkDiffusions.md)
  : Network diffusion class
- [`calculate_diffusion_auc()`](https://gregorlueg.github.io/bixverse/reference/calculate_diffusion_auc.md)
  : Calculate the AUROC for a diffusion score
- [`community_detection()`](https://gregorlueg.github.io/bixverse/reference/community_detection.md)
  : Identify privileged communities based on a given diffusion vector
- [`constrained_page_rank()`](https://gregorlueg.github.io/bixverse/reference/constrained_page_rank.md)
  : Constrained personalised page rank
- [`constrained_page_rank_ls()`](https://gregorlueg.github.io/bixverse/reference/constrained_page_rank_ls.md)
  : Constrained personalised page rank over a list
- [`diffuse_seed_nodes()`](https://gregorlueg.github.io/bixverse/reference/diffuse_seed_nodes.md)
  : Diffuse seed genes over a network
- [`permute_seed_nodes()`](https://gregorlueg.github.io/bixverse/reference/permute_seed_nodes.md)
  : Generate permuation scores for the diffusion
- [`get_diffusion_perms()`](https://gregorlueg.github.io/bixverse/reference/get_diffusion_perms.md)
  : Get the diffusion permutations
- [`generate_personalisation_vec()`](https://gregorlueg.github.io/bixverse/reference/generate_personalisation_vec.md)
  : Helper function to create personalisation vectors
- [`tied_diffusion()`](https://gregorlueg.github.io/bixverse/reference/tied_diffusion.md)
  : Diffuse seed genes in a tied manner over a network
- [`RbhGraph()`](https://gregorlueg.github.io/bixverse/reference/RbhGraph.md)
  : Reciprocal best hit graph
- [`find_rbh_communities()`](https://gregorlueg.github.io/bixverse/reference/find_rbh_communities.md)
  : Find RBH communities
- [`generate_rbh_graph()`](https://gregorlueg.github.io/bixverse/reference/generate_rbh_graph.md)
  : Generate an RBH graph.
- [`get_diffusion_vector()`](https://gregorlueg.github.io/bixverse/reference/get_diffusion_vector.md)
  : Get the diffusion vector
- [`get_rbh_res()`](https://gregorlueg.github.io/bixverse/reference/get_rbh_res.md)
  : Get the RBH results
- [`SimilarityNetworkFusion()`](https://gregorlueg.github.io/bixverse/reference/SimilarityNetworkFusion.md)
  : Similarity network fusion
- [`add_snf_data_modality()`](https://gregorlueg.github.io/bixverse/reference/add_snf_data_modality.md)
  : Add a data modality for SNF generation
- [`get_snf_params()`](https://gregorlueg.github.io/bixverse/reference/get_snf_params.md)
  : Get the SNF params
- [`get_snf_final_mat()`](https://gregorlueg.github.io/bixverse/reference/get_snf_final_mat.md)
  : Get the final SNF matrix
- [`get_snf_adjcacency_mat()`](https://gregorlueg.github.io/bixverse/reference/get_snf_adjcacency_mat.md)
  : Get an individual affinity matrix
- [`params_graph_resolution()`](https://gregorlueg.github.io/bixverse/reference/params_graph_resolution.md)
  : Wrapper function to generate resolution parameters for Leiden or
  Louvain clustering.
- [`run_snf()`](https://gregorlueg.github.io/bixverse/reference/run_snf.md)
  : Run the SNF algorithm
- [`params_community_detection()`](https://gregorlueg.github.io/bixverse/reference/params_community_detection.md)
  : Wrapper function to generate community detection parameters
- [`params_snf()`](https://gregorlueg.github.io/bixverse/reference/params_snf.md)
  : Wrapper function to generate SNF parameters

## Single cell class and getters

THE single cell class with a large number of getters.

- [`SingleCells()`](https://gregorlueg.github.io/bixverse/reference/SingleCells.md)
  : bixverse SingleCells class
- [`SingleCellCountData`](https://gregorlueg.github.io/bixverse/reference/SingleCellCountData.md)
  [`$.SingleCellCountData`](https://gregorlueg.github.io/bixverse/reference/SingleCellCountData.md)
  **\[experimental\]** : Single cell count data handler
- [`add_sc_new_obs()`](https://gregorlueg.github.io/bixverse/reference/add_sc_new_obs.md)
  : Add an obs table derived from a method to the SingleCells.
- [`get_sc_obs()`](https://gregorlueg.github.io/bixverse/reference/get_sc_obs.md)
  : Getter the obs table
- [`get_sc_var()`](https://gregorlueg.github.io/bixverse/reference/get_sc_var.md)
  : Getter the var table
- [`get_sc_counts()`](https://gregorlueg.github.io/bixverse/reference/get_sc_counts.md)
  : Getter the counts
- [`get_available_embeddings()`](https://gregorlueg.github.io/bixverse/reference/get_available_embeddings.md)
  : Get the available embeddings
- [`get_cell_indices()`](https://gregorlueg.github.io/bixverse/reference/get_cell_indices.md)
  : Get the index position for a gene
- [`get_cell_names()`](https://gregorlueg.github.io/bixverse/reference/get_cell_names.md)
  : Get the cell names
- [`get_cells_to_keep()`](https://gregorlueg.github.io/bixverse/reference/get_cells_to_keep.md)
  : Get the cells to keep
- [`get_embedding()`](https://gregorlueg.github.io/bixverse/reference/get_embedding.md)
  : Get the embedding
- [`get_gene_indices()`](https://gregorlueg.github.io/bixverse/reference/get_gene_indices.md)
  : Get the index position for a gene
- [`get_gene_names()`](https://gregorlueg.github.io/bixverse/reference/get_gene_names.md)
  : Get the gene names
- [`get_cell_info()`](https://gregorlueg.github.io/bixverse/reference/get_cell_info.md)
  : Get the cell idx (R-based) and cell names
- [`get_hvg()`](https://gregorlueg.github.io/bixverse/reference/get_hvg.md)
  : Get the HVG
- [`get_knn_mat()`](https://gregorlueg.github.io/bixverse/reference/get_knn_mat.md)
  : Get the KNN matrix
- [`get_knn_obj()`](https://gregorlueg.github.io/bixverse/reference/get_knn_obj.md)
  : Get the KNN object
- [`get_pca_singular_val()`](https://gregorlueg.github.io/bixverse/reference/get_pca_singular_val.md)
  : Get the PCA singular values
- [`get_pca_loadings()`](https://gregorlueg.github.io/bixverse/reference/get_pca_loadings.md)
  : Get the PCA loadings
- [`get_pca_factors()`](https://gregorlueg.github.io/bixverse/reference/get_pca_factors.md)
  : Get the PCA factors
- [`get_snn_graph()`](https://gregorlueg.github.io/bixverse/reference/get_snn_graph.md)
  : Get the sNN graph
- [`get_gene_names_from_idx()`](https://gregorlueg.github.io/bixverse/reference/get_gene_names_from_idx.md)
  : Get the gene names based on the gene idx
- [`get_sc_available_features()`](https://gregorlueg.github.io/bixverse/reference/get_sc_available_features.md)
  : Returns the available features for single cell applications
- [`setnames_sc()`](https://gregorlueg.github.io/bixverse/reference/setnames_sc.md)
  : Rename columns in the obs or var table
- [`set_sc_new_obs_col()`](https://gregorlueg.github.io/bixverse/reference/set_sc_new_obs_col.md)
  : Add a new column to the obs table
- [`set_sc_new_obs_col_multiple()`](https://gregorlueg.github.io/bixverse/reference/set_sc_new_obs_col_multiple.md)
  : Add multiple new columns to the obs table
- [`set_sc_new_var_cols()`](https://gregorlueg.github.io/bixverse/reference/set_sc_new_var_cols.md)
  : Add a new column to the var table
- [`drop_cols_sc()`](https://gregorlueg.github.io/bixverse/reference/drop_cols_sc.md)
  : Drop columns from the obs or var table

## Single cell class for multi-modal data and getters

This version allows you to also work ADT counts

- [`SingleCellsMultiModal()`](https://gregorlueg.github.io/bixverse/reference/SingleCellsMultiModal.md)
  : bixverse SingleCells (multi modal) class

- [`new_adt_counts_clr()`](https://gregorlueg.github.io/bixverse/reference/new_adt_counts_clr.md)
  :

  Generates a new `ADTCounts` class

- [`new_adt_counts_dsb()`](https://gregorlueg.github.io/bixverse/reference/new_adt_counts_dsb.md)
  :

  Generates a new `ADTCounts` class via DSB normalisation

- [`add_adt_counts_sc()`](https://gregorlueg.github.io/bixverse/reference/add_adt_counts_sc.md)
  :

  Add ADT counts to `SingleCellsMultiModal`

- [`detect_adt_isotypes()`](https://gregorlueg.github.io/bixverse/reference/detect_adt_isotypes.md)
  : Detect likely isotype-control features by name pattern

- [`get_adt_feature_info()`](https://gregorlueg.github.io/bixverse/reference/get_adt_feature_info.md)
  : Get the ADT feature info

- [`get_adt_names()`](https://gregorlueg.github.io/bixverse/reference/get_adt_names.md)
  : Get the ADT feature names

- [`get_adt_sample_info()`](https://gregorlueg.github.io/bixverse/reference/get_adt_sample_info.md)
  : Get the ADT sample info

- [`params_sc_dsb()`](https://gregorlueg.github.io/bixverse/reference/params_sc_dsb.md)
  : Default parameters for DSB ADT normalisation

- [`read_multi_tenx_h5_adt()`](https://gregorlueg.github.io/bixverse/reference/read_multi_tenx_h5_adt.md)
  : Read in 10x h5 ADT data from multiple files

- [`read_tenx_h5_adt()`](https://gregorlueg.github.io/bixverse/reference/read_tenx_h5_adt.md)
  : Read in 10x h5 ADT data

- [`remove_adt_isotypes()`](https://gregorlueg.github.io/bixverse/reference/remove_adt_isotypes.md)
  : Return the ADT feature names removing the isotypes

## Meta cell-related generation, classes and methods

Generating meta cells.

- [`MetaCells()`](https://gregorlueg.github.io/bixverse/reference/MetaCells.md)
  : bixverse meta cell class

- [`calc_diffusion_coordinates()`](https://gregorlueg.github.io/bixverse/reference/calc_diffusion_coordinates.md)
  : Calculate diffusion coordinates

- [`calc_manifold_metrics()`](https://gregorlueg.github.io/bixverse/reference/calc_manifold_metrics.md)
  : Calculate manifold metrics

- [`calc_meta_cell_purity()`](https://gregorlueg.github.io/bixverse/reference/calc_meta_cell_purity.md)
  : Calculate meta cell purity

- [`generate_supercells_sc()`](https://gregorlueg.github.io/bixverse/reference/generate_supercells_sc.md)
  :

  Generate SuperCells and return a `MetaCells` object

- [`generate_bt_meta_cells_sc()`](https://gregorlueg.github.io/bixverse/reference/generate_bt_meta_cells_sc.md)
  :

  Generate meta cells based on hdWGCNA and return a `MetaCells` object

- [`generate_seacells_sc()`](https://gregorlueg.github.io/bixverse/reference/generate_seacells_sc.md)
  :

  Generate meta cells based on SEACells and return a `MetaCells` object

- [`params_sc_supercell()`](https://gregorlueg.github.io/bixverse/reference/params_sc_supercell.md)
  : Wrapper function for parameters for SuperCell generation

- [`params_sc_bt_metacells()`](https://gregorlueg.github.io/bixverse/reference/params_sc_bt_metacells.md)
  : Wrapper function for parameters for bootstrapped meta cell
  generation

- [`params_sc_seacells()`](https://gregorlueg.github.io/bixverse/reference/params_sc_seacells.md)
  : Wrapper function for the SEACells parameters

## Single cell i/o

I/O functions for single cell. Load in h5ad, mtx, Seurat or R data into
Rust and the DuckDB supporting the metadata.

- [`get_cell_ranger_params()`](https://gregorlueg.github.io/bixverse/reference/get_cell_ranger_params.md)
  : Helper to generate cell ranger input parameters

- [`get_h5ad_dimensions()`](https://gregorlueg.github.io/bixverse/reference/get_h5ad_dimensions.md)
  : Helper function to get the dimensions and storage format

- [`prescan_h5ad_files()`](https://gregorlueg.github.io/bixverse/reference/prescan_h5ad_files.md)
  : Pre-scan multiple h5ad files for multi-sample loading

- [`prescan_mtx_dirs()`](https://gregorlueg.github.io/bixverse/reference/prescan_mtx_dirs.md)
  : Prescan multiple mtx directories for a multi-load

- [`prescan_tenx_h5_files()`](https://gregorlueg.github.io/bixverse/reference/prescan_tenx_h5_files.md)
  : Pre-scan multiple 10x CellRanger h5 files for multi-sample loading

- [`read_tenx_h5_metadata()`](https://gregorlueg.github.io/bixverse/reference/read_tenx_h5_metadata.md)
  : Read barcode and feature tables and metadata from a 10x h5 file

- [`load_existing()`](https://gregorlueg.github.io/bixverse/reference/load_existing.md)
  : Load an existing SingleCells from disk

- [`load_h5ad()`](https://gregorlueg.github.io/bixverse/reference/load_h5ad.md)
  :

  Load in h5ad to `SingleCells`

- [`load_h5ad_norm()`](https://gregorlueg.github.io/bixverse/reference/load_h5ad_norm.md)
  :

  Load in h5ad with normalised counts to `SingleCells`

- [`load_mtx()`](https://gregorlueg.github.io/bixverse/reference/load_mtx.md)
  :

  Load in mtx/plain text files to `SingleCells`

- [`load_multi_mtx()`](https://gregorlueg.github.io/bixverse/reference/load_multi_mtx.md)
  :

  Load multiple mtx directories into a single `SingleCells`

- [`load_multi_h5ad()`](https://gregorlueg.github.io/bixverse/reference/load_multi_h5ad.md)
  :

  Load multiple h5ad files into a single `SingleCells`

- [`stream_h5ad()`](https://gregorlueg.github.io/bixverse/reference/stream_h5ad.md)
  :

  Stream in h5ad to `SingleCells` (alias)

- [`load_r_data()`](https://gregorlueg.github.io/bixverse/reference/load_r_data.md)
  : Load in data directly from R objects.

- [`load_seurat()`](https://gregorlueg.github.io/bixverse/reference/load_seurat.md)
  :

  Load in Seurat to `SingleCells`

- [`load_tenx_h5()`](https://gregorlueg.github.io/bixverse/reference/load_tenx_h5.md)
  :

  Load in a 10x CellRanger h5 file to `SingleCells`

- [`load_multi_tenx_h5()`](https://gregorlueg.github.io/bixverse/reference/load_multi_tenx_h5.md)
  :

  Load multiple 10x CellRanger h5 files into a single `SingleCells`

- [`read_h5ad_metadata()`](https://gregorlueg.github.io/bixverse/reference/read_h5ad_metadata.md)
  : Read obs and var tables and metadata from an h5ad file

- [`read_h5ad_x_summary()`](https://gregorlueg.github.io/bixverse/reference/read_h5ad_x_summary.md)
  : Read summary statistics from the X slot of an h5ad file

- [`save_sc_exp_to_disk()`](https://gregorlueg.github.io/bixverse/reference/save_sc_exp_to_disk.md)
  : Save memory-bound data to disk

- [`sc_old_file_conversion()`](https://gregorlueg.github.io/bixverse/reference/sc_old_file_conversion.md)
  : Convert legacy v2 single-cell data files to v3 format

- [`merge_sc_experiments()`](https://gregorlueg.github.io/bixverse/reference/merge_sc_experiments.md)
  :

  Merge multiple `SingleCells` experiments into one

- [`params_sc_min_quality()`](https://gregorlueg.github.io/bixverse/reference/params_sc_min_quality.md)
  : Wrapper function to generate QC metric params for single cell

- [`params_sc_mtx_io()`](https://gregorlueg.github.io/bixverse/reference/params_sc_mtx_io.md)
  : Wrapper function to provide data for mtx-based loading

## Single cell processing

Helpers to process single cell data. Doublet detection, proportions of
gene sets, HVG (batch-aware), PCA and batch corrections.

- [`scrublet_sc()`](https://gregorlueg.github.io/bixverse/reference/scrublet_sc.md)
  : Doublet detection with Scrublet
- [`call_doublets_manual()`](https://gregorlueg.github.io/bixverse/reference/call_doublets_manual.md)
  : Manually readjust Scrublet doublet call thresholds
- [`doublet_detection_boost_sc()`](https://gregorlueg.github.io/bixverse/reference/doublet_detection_boost_sc.md)
  : Doublet detection with boosted doublet classification
- [`scdblfinder_sc()`](https://gregorlueg.github.io/bixverse/reference/scdblfinder_sc.md)
  : Run scDblFinder doublet detection on a SingleCells object
- [`gene_set_proportions_sc()`](https://gregorlueg.github.io/bixverse/reference/gene_set_proportions_sc.md)
  : Calculate the proportions of reads for specific gene sets
- [`per_cell_qc_outlier()`](https://gregorlueg.github.io/bixverse/reference/per_cell_qc_outlier.md)
  : Use MAD outlier detection on per-cell QC metrics
- [`run_cell_qc()`](https://gregorlueg.github.io/bixverse/reference/run_cell_qc.md)
  : Run MAD outlier detection on per-cell QC metrics
- [`find_hvg_sc()`](https://gregorlueg.github.io/bixverse/reference/find_hvg_sc.md)
  : Identify HVGs
- [`find_hvg_batch_aware_sc()`](https://gregorlueg.github.io/bixverse/reference/find_hvg_batch_aware_sc.md)
  : Identify HVGs (batch aware)
- [`get_hvg_data_sc()`](https://gregorlueg.github.io/bixverse/reference/get_hvg_data_sc.md)
  : Identify HVGs without mutating object state
- [`calculate_pca_sc()`](https://gregorlueg.github.io/bixverse/reference/calculate_pca_sc.md)
  : Run PCA for single cell
- [`generate_sc_knn()`](https://gregorlueg.github.io/bixverse/reference/generate_sc_knn.md)
  : Generate a new SingleCellNearestNeighbour from data
- [`find_neighbours_sc()`](https://gregorlueg.github.io/bixverse/reference/find_neighbours_sc.md)
  : Find the neighbours for single cell.
- [`top_genes_perc_sc()`](https://gregorlueg.github.io/bixverse/reference/top_genes_perc_sc.md)
  : Calculate the proportions of reads for the Top N genes
- [`fast_mnn_sc()`](https://gregorlueg.github.io/bixverse/reference/fast_mnn_sc.md)
  : Run fastMNN
- [`harmony_sc()`](https://gregorlueg.github.io/bixverse/reference/harmony_sc.md)
  : Run Harmony
- [`harmony_v2_sc()`](https://gregorlueg.github.io/bixverse/reference/harmony_v2_sc.md)
  : Run Harmony v2
- [`bbknn_sc()`](https://gregorlueg.github.io/bixverse/reference/bbknn_sc.md)
  : Run BBKNN
- [`calculate_kbet_sc()`](https://gregorlueg.github.io/bixverse/reference/calculate_kbet_sc.md)
  : Calculate kBET scores
- [`calculate_batch_asw_sc()`](https://gregorlueg.github.io/bixverse/reference/calculate_batch_asw_sc.md)
  : Calculate batch average silhouette width
- [`calculate_batch_lisi_sc()`](https://gregorlueg.github.io/bixverse/reference/calculate_batch_lisi_sc.md)
  : Calculate batch LISI scores
- [`params_norm_doublets_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_norm_doublets_defaults.md)
  : Helper function to generate normalisation defaults for doublet
  detection.
- [`params_boost()`](https://gregorlueg.github.io/bixverse/reference/params_boost.md)
  : Wrapper function for Boost parameters
- [`params_sc_hvg()`](https://gregorlueg.github.io/bixverse/reference/params_sc_hvg.md)
  : Wrapper function for HVG detection parameters.
- [`params_sc_pca()`](https://gregorlueg.github.io/bixverse/reference/params_sc_pca.md)
  : Wrapper for PCA specifically designed for single cells
- [`params_scrublet()`](https://gregorlueg.github.io/bixverse/reference/params_scrublet.md)
  : Wrapper function for Scrublet doublet detection parameters
- [`params_sc_bbknn()`](https://gregorlueg.github.io/bixverse/reference/params_sc_bbknn.md)
  : Wrapper function for the BBKNN parameters
- [`params_sc_fast_cluster()`](https://gregorlueg.github.io/bixverse/reference/params_sc_fast_cluster.md)
  : Fast single cell clustering parameters
- [`params_sc_fastmnn()`](https://gregorlueg.github.io/bixverse/reference/params_sc_fastmnn.md)
  : Wrapper function for the fastMNN parameters
- [`params_sc_harmony()`](https://gregorlueg.github.io/bixverse/reference/params_sc_harmony.md)
  : Default parameters for Harmony batch correction
- [`params_sc_harmony_v2()`](https://gregorlueg.github.io/bixverse/reference/params_sc_harmony_v2.md)
  : Default parameters for Harmony v2 batch correction
- [`params_sc_neighbours()`](https://gregorlueg.github.io/bixverse/reference/params_sc_neighbours.md)
  : Wrapper function for parameters for neighbour identification in
  single cell
- [`params_scdblfinder()`](https://gregorlueg.github.io/bixverse/reference/params_scdblfinder.md)
  : Wrapper function for scDblFinder doublet detection parameters
- [`params_hvg_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_hvg_defaults.md)
  : Helper function to generate HVG defaults
- [`params_pca_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_pca_defaults.md)
  : Helper function to generate default parameters for PCA
- [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
  : Helper function to generate kNN defaults
- [`params_sc_knn()`](https://gregorlueg.github.io/bixverse/reference/params_sc_knn.md)
  : Parameters for single cell kNN searches
- [`params_kmeans_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_kmeans_defaults.md)
  : K-mean parameter defaults.
- [`params_fast_cluster_default()`](https://gregorlueg.github.io/bixverse/reference/params_fast_cluster_default.md)
  : Helper function to generate default parameters for the fast
  clustering for the doublet detection methods

## Single cell analysis methods

A large number of different methods to extract insights from your single
cell experiment. Gene set scoring, DGEs, kNN generations, pseudo-bulk
count extraction, miloR, Hotspot, VISION and SCENIC.

- [`aucell_sc()`](https://gregorlueg.github.io/bixverse/reference/aucell_sc.md)
  : Calculate AUC scores (akin to AUCell)

- [`module_scores_sc()`](https://gregorlueg.github.io/bixverse/reference/module_scores_sc.md)
  : Calculate module activity scores

- [`fast_cluster_sc()`](https://gregorlueg.github.io/bixverse/reference/fast_cluster_sc.md)
  : Run fast Louvain clustering on a SingleCells object

- [`find_clusters_sc()`](https://gregorlueg.github.io/bixverse/reference/find_clusters_sc.md)
  : Graph-based clustering of cells on the sNN graph

- [`find_markers_sc()`](https://gregorlueg.github.io/bixverse/reference/find_markers_sc.md)
  : Calculate DGE between two cell groups

- [`find_all_markers_sc()`](https://gregorlueg.github.io/bixverse/reference/find_all_markers_sc.md)
  : Find all markers

- [`get_pseudobulked_sc()`](https://gregorlueg.github.io/bixverse/reference/get_pseudobulked_sc.md)
  : Generate pseudo-bulked matrices

- [`generate_knn_sc()`](https://gregorlueg.github.io/bixverse/reference/generate_knn_sc.md)
  :

  Generate a `SingleCellNearestNeighbour` from a single cell class

- [`get_differential_abundance_res()`](https://gregorlueg.github.io/bixverse/reference/get_differential_abundance_res.md)
  : Get the differential abundance results

- [`hotspot_autocor_sc()`](https://gregorlueg.github.io/bixverse/reference/hotspot_autocor_sc.md)
  : Calculate the local auto-correlation of a gene

- [`hotspot_gene_cor_sc()`](https://gregorlueg.github.io/bixverse/reference/hotspot_gene_cor_sc.md)
  : Calculate the local pairwise gene-gene correlation

- [`generate_hotspot_membership()`](https://gregorlueg.github.io/bixverse/reference/generate_hotspot_membership.md)
  : Identify hotspot gene clusters

- [`get_hotspot_membership()`](https://gregorlueg.github.io/bixverse/reference/get_hotspot_membership.md)
  : Get the hotspot gene membership table

- [`get_miloR_abundances_sc()`](https://gregorlueg.github.io/bixverse/reference/get_miloR_abundances_sc.md)
  : Generate an miloR abundance object for differential abundance
  testing

- [`meld_sc()`](https://gregorlueg.github.io/bixverse/reference/meld_sc.md)
  : Run MELD signal smoothing for differential abundance estimation

- [`get_index_cells()`](https://gregorlueg.github.io/bixverse/reference/get_index_cells.md)
  : Get the index cells

- [`add_nhoods_info()`](https://gregorlueg.github.io/bixverse/reference/add_nhoods_info.md)
  : Add neighbourhood info on majority cell type

- [`test_nhoods()`](https://gregorlueg.github.io/bixverse/reference/test_nhoods.md)
  : Test neighbourhoods for differential abundance

- [`vision_sc()`](https://gregorlueg.github.io/bixverse/reference/vision_sc.md)
  : Calculate VISION scores

- [`vision_w_autocor_sc()`](https://gregorlueg.github.io/bixverse/reference/vision_w_autocor_sc.md)
  : Calculate VISION scores (with auto-correlation scores)

- [`identify_tf_to_genes()`](https://gregorlueg.github.io/bixverse/reference/identify_tf_to_genes.md)
  : Identify the TF to gene regulation

- [`scenic_gene_filter_sc()`](https://gregorlueg.github.io/bixverse/reference/scenic_gene_filter_sc.md)
  : Filter genes for SCENIC GRN inference

- [`scenic_grn_sc()`](https://gregorlueg.github.io/bixverse/reference/scenic_grn_sc.md)
  : Run SCENIC GRN inference

- [`get_cistarget_res()`](https://gregorlueg.github.io/bixverse/reference/get_cistarget_res.md)
  : Extract the TF to gene data from the ScenicGrn object

- [`get_tf_to_gene()`](https://gregorlueg.github.io/bixverse/reference/get_tf_to_gene.md)
  : Extract the TF to gene data from the ScenicGrn object

- [`tf_to_genes_correlations()`](https://gregorlueg.github.io/bixverse/reference/tf_to_genes_correlations.md)
  : Generate TF to gene correlations

- [`tf_to_genes_motif_enrichment()`](https://gregorlueg.github.io/bixverse/reference/tf_to_genes_motif_enrichment.md)
  : Run the SCENIC motif enrichment

- [`nmf_sc()`](https://gregorlueg.github.io/bixverse/reference/nmf_sc.md)
  : Run single-run NMF on single cell or meta cell data

- [`stabilised_nmf_sc()`](https://gregorlueg.github.io/bixverse/reference/stabilised_nmf_sc.md)
  : Run stabilised (multi-run) NMF on single cell or meta cell data

- [`get_best_run()`](https://gregorlueg.github.io/bixverse/reference/get_best_run.md)
  : Get the best run from a stabilised NMF result

- [`get_w()`](https://gregorlueg.github.io/bixverse/reference/get_w.md)
  : Get the W (gene loadings) matrix

- [`get_h()`](https://gregorlueg.github.io/bixverse/reference/get_h.md)
  : Get the H (cell activations) matrix

- [`params_sc_hotspot()`](https://gregorlueg.github.io/bixverse/reference/params_sc_hotspot.md)
  : Wrapper function for parameters for HotSpot

- [`params_sc_miloR()`](https://gregorlueg.github.io/bixverse/reference/params_sc_miloR.md)
  : Wrapper function for parameters for MiloR

- [`params_sc_vision()`](https://gregorlueg.github.io/bixverse/reference/params_sc_vision.md)
  : Wrapper function for parameters for VISION with auto-correlation

- [`params_scenic()`](https://gregorlueg.github.io/bixverse/reference/params_scenic.md)
  : Constructor for SCENIC parameters

- [`params_scenic_extra_trees_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_scenic_extra_trees_defaults.md)
  : Default parameters for the SCENIC ExtraTrees regression learner

- [`params_scenic_gradient_boosting_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_scenic_gradient_boosting_defaults.md)
  : Default parameters for the SCENIC GradientBoosting (GRNBoost2)
  regression learner

- [`params_scenic_random_forest_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_scenic_random_forest_defaults.md)
  : Default parameters for the SCENIC RandomForest regression learner

- [`params_meld()`](https://gregorlueg.github.io/bixverse/reference/params_meld.md)
  : Constructor for MELD parameters

- [`params_nmf_hals()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_hals.md)
  : Wrapper function for NMF (HALS) parameters

## Single cell multi-modal analysis methods

Methods to analyse multi-modal single cell data

- [`calculate_pca_adt_sc()`](https://gregorlueg.github.io/bixverse/reference/calculate_pca_adt_sc.md)
  : Calculate the PCA on top of the normalised ADT counts
- [`generate_wnn_graph_sc()`](https://gregorlueg.github.io/bixverse/reference/generate_wnn_graph_sc.md)
  : Generate a weighted nearest neighbour (WNN) graph
- [`params_sc_wnn()`](https://gregorlueg.github.io/bixverse/reference/params_sc_wnn.md)
  : Wrapper function for WNN parameters

## Single-cell related classes and methods

Additional helpers for specific small sub classes used in single cell.

- [`calc_knn_metrics()`](https://gregorlueg.github.io/bixverse/reference/calc_knn_metrics.md)
  : Calculate recall at k and distance ratio
- [`get_centroids()`](https://gregorlueg.github.io/bixverse/reference/get_centroids.md)
  : Get k-means centroids from a fast cluster result
- [`get_feature_mat()`](https://gregorlueg.github.io/bixverse/reference/get_feature_mat.md)
  : Get the feature matrix used for the classifier
- [`get_kmeans_clusters()`](https://gregorlueg.github.io/bixverse/reference/get_kmeans_clusters.md)
  : Get k-means cluster assignments from a fast cluster result
- [`get_knn_dist()`](https://gregorlueg.github.io/bixverse/reference/get_knn_dist.md)
  : Get the KNN distance
- [`get_data()`](https://gregorlueg.github.io/bixverse/reference/get_data.md)
  : Get the ready obs data from various sub method
- [`get_scores()`](https://gregorlueg.github.io/bixverse/reference/get_scores.md)
  : Get scores
- [`new_sc_knn()`](https://gregorlueg.github.io/bixverse/reference/new_sc_knn.md)
  : Helper function to generate kNN data with distances

## Reference mapping and cell type annotations in single cell

Helpers to do reference mapping and cell type identification in single
cell

- [`SymphonyReference()`](https://gregorlueg.github.io/bixverse/reference/SymphonyReference.md)
  : bixverse SymphonyReference class
- [`add_symphony_labels()`](https://gregorlueg.github.io/bixverse/reference/add_symphony_labels.md)
  : Add labels to a Symphony reference post-hoc
- [`build_symphony_ref()`](https://gregorlueg.github.io/bixverse/reference/build_symphony_ref.md)
  : Build a Symphony reference from a SingleCells object
- [`transfer_labels_symphony()`](https://gregorlueg.github.io/bixverse/reference/transfer_labels_symphony.md)
  : Transfer labels from a Symphony reference to a query via kNN
  majority vote
- [`get_symphony_hvg_names()`](https://gregorlueg.github.io/bixverse/reference/get_symphony_hvg_names.md)
  : Getter for the HVG gene names of a Symphony reference
- [`get_symphony_labels()`](https://gregorlueg.github.io/bixverse/reference/get_symphony_labels.md)
  : Getter for the stored labels of a Symphony reference
- [`get_symphony_loadings()`](https://gregorlueg.github.io/bixverse/reference/get_symphony_loadings.md)
  : Getter for the PCA loadings of a Symphony reference
- [`get_symphony_z_corr()`](https://gregorlueg.github.io/bixverse/reference/get_symphony_z_corr.md)
  : Getter for the corrected embedding of a Symphony reference
- [`map_symphony_query()`](https://gregorlueg.github.io/bixverse/reference/map_symphony_query.md)
  : Map a SingleCells query onto a Symphony reference
- [`prepare_cell_markers()`](https://gregorlueg.github.io/bixverse/reference/prepare_cell_markers.md)
  : Helper function to prepare cell markers
- [`calc_sc_type_scores()`](https://gregorlueg.github.io/bixverse/reference/calc_sc_type_scores.md)
  : Calculate ScType scores per cell
- [`score_clusters()`](https://gregorlueg.github.io/bixverse/reference/score_clusters.md)
  : Score clusters based on ScType
- [`params_symphony_map()`](https://gregorlueg.github.io/bixverse/reference/params_symphony_map.md)
  : Default parameters for Symphony query mapping

## Single cell ligand receptor

Classes, functions and generics/methods for ligand receptor analysis for
single cell.

- [`compute_expression_info_sc()`](https://gregorlueg.github.io/bixverse/reference/compute_expression_info_sc.md)
  : Compute per-cluster mean expression and expressing fraction for a
  gene set
- [`generate_ligand_target_influence()`](https://gregorlueg.github.io/bixverse/reference/generate_ligand_target_influence.md)
  : Generate the ligand to target influence matrix
- [`get_influence()`](https://gregorlueg.github.io/bixverse/reference/get_influence.md)
  : Get the ligand-target influence matrix
- [`ligand_activity_scores()`](https://gregorlueg.github.io/bixverse/reference/ligand_activity_scores.md)
  : Compute ligand activity scores against gene sets
- [`params_ligand_target()`](https://gregorlueg.github.io/bixverse/reference/params_ligand_target.md)
  : Parameters for ligand to target influence computation
- [`prioritise_interactions()`](https://gregorlueg.github.io/bixverse/reference/prioritise_interactions.md)
  : Prioritise sender-ligand-receiver-receptor interactions

## Single cell plotting stuff

Various helpers that generate data for plotting single cell, such as 2D
embeddings, or extract specific columns from the metadata or genes (and
their summaries) from the binary storage files.

- [`sc_knn_to_nearest_neighbours()`](https://gregorlueg.github.io/bixverse/reference/sc_knn_to_nearest_neighbours.md)
  : Convert SingleCellNearestNeighbour to manifoldsR NearestNeighbours
- [`umap_sc()`](https://gregorlueg.github.io/bixverse/reference/umap_sc.md)
  : Run UMAP on a SingleCells/MetaCells object
- [`tsne_sc()`](https://gregorlueg.github.io/bixverse/reference/tsne_sc.md)
  : Run t-SNE on a SingleCells/MetaCells object
- [`phate_sc()`](https://gregorlueg.github.io/bixverse/reference/phate_sc.md)
  : Run PHATE on a SingleCells/MetaCells object
- [`extract_dot_plot_data()`](https://gregorlueg.github.io/bixverse/reference/extract_dot_plot_data.md)
  : Extract grouped gene statistics for dot plots
- [`extract_gene_expression()`](https://gregorlueg.github.io/bixverse/reference/extract_gene_expression.md)
  : Extract normalised gene expression for plotting
- [`extract_embedding_data()`](https://gregorlueg.github.io/bixverse/reference/extract_embedding_data.md)
  : Extract embedding coordinates for plotting
- [`extract_feature_pair()`](https://gregorlueg.github.io/bixverse/reference/extract_feature_pair.md)
  : Extract a pair of features for scatter / hex plots
- [`extract_feature_plot_data()`](https://gregorlueg.github.io/bixverse/reference/extract_feature_plot_data.md)
  : Extract per-cell expression mapped onto an embedding
- [`extract_gene_violin_data()`](https://gregorlueg.github.io/bixverse/reference/extract_gene_violin_data.md)
  : Extract per-cell expression grouped for violin plots

## Statistical functions

Any types of functions that help with statistics

- [`calculate_effect_size()`](https://gregorlueg.github.io/bixverse/reference/calculate_effect_size.md)
  : Calculate the Hedge G effect between two matrices
- [`calculate_tom()`](https://gregorlueg.github.io/bixverse/reference/calculate_tom.md)
  : Calculate the TOM from a correlation matrix
- [`calculate_tom_from_exp()`](https://gregorlueg.github.io/bixverse/reference/calculate_tom_from_exp.md)
  : Calculate the TOM from an expression matrix
- [`f1_score_confusion_mat()`](https://gregorlueg.github.io/bixverse/reference/f1_score_confusion_mat.md)
  : F1 scores on top of a confusion matrix
- [`fast_ica_rust()`](https://gregorlueg.github.io/bixverse/reference/fast_ica_rust.md)
  : Fast ICA via Rust
- [`fast_ica_rust_helper()`](https://gregorlueg.github.io/bixverse/reference/fast_ica_rust_helper.md)
  : Fast ICA via Rust from processed data
- [`get_inflection_point()`](https://gregorlueg.github.io/bixverse/reference/get_inflection_point.md)
  : Identify the inflection point for elbow-like data
- [`ot_harmonic_score()`](https://gregorlueg.github.io/bixverse/reference/ot_harmonic_score.md)
  : Calculates a harmonic sum normalised between 0 to 1.
- [`robust_scale()`](https://gregorlueg.github.io/bixverse/reference/robust_scale.md)
  : Robust scaler.

## Plotting helpers

Some core plotting helpers in the package (usually for QC). The ones to
plot downstream results can be found in bixverse.plots.

- [`plot_boxplot_normalization()`](https://gregorlueg.github.io/bixverse/reference/plot_boxplot_normalization.md)
  : Helper plot function for boxplot of normalized data
- [`plot_epsilon_res()`](https://gregorlueg.github.io/bixverse/reference/plot_epsilon_res.md)
  : Plot the epsilon vs. power law goodness of fit result
- [`plot_hvgs()`](https://gregorlueg.github.io/bixverse/reference/plot_hvgs.md)
  : Plot the highly variable genes
- [`plot_ica_ncomp_params()`](https://gregorlueg.github.io/bixverse/reference/plot_ica_ncomp_params.md)
  : Plot various parameters with no comp
- [`plot_ica_stability_individual()`](https://gregorlueg.github.io/bixverse/reference/plot_ica_stability_individual.md)
  : Plot the stability of the ICA components
- [`plot_optimal_cuts()`](https://gregorlueg.github.io/bixverse/reference/plot_optimal_cuts.md)
  : Plot the k cuts vs median R2
- [`plot_pca()`](https://gregorlueg.github.io/bixverse/reference/plot_pca.md)
  : Helper plot function for pca with contrasts
- [`plot_pca_res()`](https://gregorlueg.github.io/bixverse/reference/plot_pca_res.md)
  : Plot the PCA data
- [`plot_preprocessing_genes()`](https://gregorlueg.github.io/bixverse/reference/plot_preprocessing_genes.md)
  : Helper plot function of distribution of genes by samples
- [`plot_preprocessing_outliers()`](https://gregorlueg.github.io/bixverse/reference/plot_preprocessing_outliers.md)
  : Helper plot function for identification of outliers
- [`plot_rbf_impact()`](https://gregorlueg.github.io/bixverse/reference/plot_rbf_impact.md)
  : Helper function to plot distance to affinity relationship
- [`plot_resolution_res()`](https://gregorlueg.github.io/bixverse/reference/plot_resolution_res.md)
  : Plot the resolution results.
- [`plot_voom_normalization()`](https://gregorlueg.github.io/bixverse/reference/plot_voom_normalization.md)
  : Helper plot function for Voom normalisation

## Data downloads and synthetic data generation

Functions and helpers to download or generate synthetic data.

- [`download_cd34_data()`](https://gregorlueg.github.io/bixverse/reference/download_cd34_data.md)
  : Download the CD34 example data from SEACells
- [`download_pbmc3k()`](https://gregorlueg.github.io/bixverse/reference/download_pbmc3k.md)
  : Download PBMC3K data from Zenodo
- [`download_demuxlet_pbmc()`](https://gregorlueg.github.io/bixverse/reference/download_demuxlet_pbmc.md)
  : Download PBMCs with demuxlet doublet information
- [`download_pbmc_batches()`](https://gregorlueg.github.io/bixverse/reference/download_pbmc_batches.md)
  : Download two different PBMC data sets for batch correction testing
- [`download_pbmc_totalseq_data()`](https://gregorlueg.github.io/bixverse/reference/download_pbmc_totalseq_data.md)
  : Download the PBMC TotalSeq data with ADT counts
- [`download_pbmc8k()`](https://gregorlueg.github.io/bixverse/reference/download_pbmc8k.md)
  : Download PBMC8K data from Zenodo
- [`calculate_sparsity_stats()`](https://gregorlueg.github.io/bixverse/reference/calculate_sparsity_stats.md)
  : Helper function to calculate the induced sparsity
- [`generate_gene_module_data()`](https://gregorlueg.github.io/bixverse/reference/generate_gene_module_data.md)
  : Generates synthetic gene module data.
- [`generate_single_cell_test_data()`](https://gregorlueg.github.io/bixverse/reference/generate_single_cell_test_data.md)
  : Single cell test data
- [`cell_cycle_genes`](https://gregorlueg.github.io/bixverse/reference/cell_cycle_genes.md)
  : Cell cycle genes
- [`write_cellranger_output()`](https://gregorlueg.github.io/bixverse/reference/write_cellranger_output.md)
  : Helper function to write data to a cell ranger like output
- [`write_h5ad_sc()`](https://gregorlueg.github.io/bixverse/reference/write_h5ad_sc.md)
  : Helper function to write data to h5ad format
- [`write_h5ad_sc_dense()`](https://gregorlueg.github.io/bixverse/reference/write_h5ad_sc_dense.md)
  : Helper function to write data to a dense h5ad file
- [`params_sc_synthetic_data()`](https://gregorlueg.github.io/bixverse/reference/params_sc_synthetic_data.md)
  : Default parameters for generation of synthetic single cell data
  (RNA)
- [`synthetic_signal_matrix()`](https://gregorlueg.github.io/bixverse/reference/synthetic_signal_matrix.md)
  : Generates a simple synthetic, pseudo gene expression matrix
- [`simulate_dropouts()`](https://gregorlueg.github.io/bixverse/reference/simulate_dropouts.md)
  : Simulate dropouts via different functions on synthetic bulk data
- [`synthetic_bulk_cor_matrix()`](https://gregorlueg.github.io/bixverse/reference/synthetic_bulk_cor_matrix.md)
  : Generates synthetic bulk RNAseq data
- [`synthetic_c_pca_data()`](https://gregorlueg.github.io/bixverse/reference/synthetic_c_pca_data.md)
  : Generates synthetic data for contrastive PCA exploration.

## Utils

All types of other random helpers without a clear pattern

- [`AnnDataParser`](https://gregorlueg.github.io/bixverse/reference/AnnDataParser.md)
  : Class for Anndata
- [`calculate_sparsity_stats()`](https://gregorlueg.github.io/bixverse/reference/calculate_sparsity_stats.md)
  : Helper function to calculate the induced sparsity
- [`find_threshold_otsu()`](https://gregorlueg.github.io/bixverse/reference/find_threshold_otsu.md)
  : Find a threshold via the Otsu method
- [`get_seurat_counts_to_list()`](https://gregorlueg.github.io/bixverse/reference/get_seurat_counts_to_list.md)
  : Transform Seurat raw counts into a List
- [`knn_graph_label_propagation()`](https://gregorlueg.github.io/bixverse/reference/knn_graph_label_propagation.md)
  : kNN-based graph label propagation
- [`params_label_propagation()`](https://gregorlueg.github.io/bixverse/reference/params_label_propagation.md)
  : Wrapper function to generate label propagation parameters
- [`upper_triangle_to_sparse()`](https://gregorlueg.github.io/bixverse/reference/upper_triangle_to_sparse.md)
  : Transform an upper triangle-stored matrix to a sparse one
- [`upper_triangular_sym_mat`](https://gregorlueg.github.io/bixverse/reference/upper_triangular_sym_mat.md)
  : Class for symmetric correlation matrices
- [`to_snake_case()`](https://gregorlueg.github.io/bixverse/reference/to_snake_case.md)
  : Helper function to transform strings to snake_case

## Rust wrappers

Everything Rusty - only use this if you know what you are doing… Maybe
useful for your own package? Use with care and read the documentation!
The ones exposed here are general enough to be useful in other packages.
There is a lot more under the hood…

- [`rs_2d_loess()`](https://gregorlueg.github.io/bixverse/reference/rs_2d_loess.md)
  **\[experimental\]** : Rust implementation of a Loess function
- [`rs_batch_lisi()`](https://gregorlueg.github.io/bixverse/reference/rs_batch_lisi.md)
  **\[experimental\]** : Calculate batch LISI scores
- [`rs_batch_silhouette_width()`](https://gregorlueg.github.io/bixverse/reference/rs_batch_silhouette_width.md)
  **\[experimental\]** : Calculate batch silhouette width from an
  embedding
- [`rs_cistarget()`](https://gregorlueg.github.io/bixverse/reference/rs_cistarget.md)
  **\[experimental\]** : Run CisTarget motif enrichment analysis
- [`rs_compare_knn()`](https://gregorlueg.github.io/bixverse/reference/rs_compare_knn.md)
  **\[experimental\]** : Helper to compare kNN graphs
- [`rs_constrained_page_rank()`](https://gregorlueg.github.io/bixverse/reference/rs_constrained_page_rank.md)
  **\[experimental\]** : Calculate a constrained page-rank score
- [`rs_constrained_page_rank_list()`](https://gregorlueg.github.io/bixverse/reference/rs_constrained_page_rank_list.md)
  **\[experimental\]** : Calculate a constrained page-rank score over a
  list.
- [`rs_contrastive_pca()`](https://gregorlueg.github.io/bixverse/reference/rs_contrastive_pca.md)
  **\[experimental\]** : Calculate the contrastive PCA
- [`rs_cor()`](https://gregorlueg.github.io/bixverse/reference/rs_cor.md)
  **\[experimental\]** : Calculate the column wise correlations.
- [`rs_cor2()`](https://gregorlueg.github.io/bixverse/reference/rs_cor2.md)
  **\[experimental\]** : Calculate the column wise correlations.
- [`rs_cor_upper_triangle()`](https://gregorlueg.github.io/bixverse/reference/rs_cor_upper_triangle.md)
  **\[experimental\]** : Calculate the column wise correlations and
  returns the upper triangle
- [`rs_cos()`](https://gregorlueg.github.io/bixverse/reference/rs_cos.md)
  **\[experimental\]** : Calculate the column wise cosine similarities
- [`rs_count_zeroes()`](https://gregorlueg.github.io/bixverse/reference/rs_count_zeroes.md)
  **\[experimental\]** : Helper to get zero stats from a given matrix
- [`rs_cov2cor()`](https://gregorlueg.github.io/bixverse/reference/rs_cov2cor.md)
  **\[experimental\]** : Calculates the correlation matrix from the
  co-variance matrix
- [`rs_covariance()`](https://gregorlueg.github.io/bixverse/reference/rs_covariance.md)
  **\[experimental\]** : Calculate the column-wise co-variance.
- [`rs_dense_to_upper_triangle()`](https://gregorlueg.github.io/bixverse/reference/rs_dense_to_upper_triangle.md)
  **\[experimental\]** : Generate a vector-based representation of the
  upper triangle of a matrix
- [`rs_differential_cor()`](https://gregorlueg.github.io/bixverse/reference/rs_differential_cor.md)
  **\[experimental\]** : Calculate the column wise differential
  correlation between two sets of data.
- [`rs_dist()`](https://gregorlueg.github.io/bixverse/reference/rs_dist.md)
  **\[experimental\]** : Calculate the pairwise column distance in a
  matrix
- [`rs_fast_auc()`](https://gregorlueg.github.io/bixverse/reference/rs_fast_auc.md)
  **\[experimental\]** : Fast AUC calculation
- [`rs_fast_ica()`](https://gregorlueg.github.io/bixverse/reference/rs_fast_ica.md)
  **\[experimental\]** : Run the Rust implementation of fast ICA.
- [`rs_fdr_adjustment()`](https://gregorlueg.github.io/bixverse/reference/rs_fdr_adjustment.md)
  **\[experimental\]** : Calculate a BH-based FDR
- [`rs_geom_elim_fgsea_simple()`](https://gregorlueg.github.io/bixverse/reference/rs_geom_elim_fgsea_simple.md)
  **\[experimental\]** : Run fgsea simple method for gene ontology with
  elimination method
- [`rs_gower_dist()`](https://gregorlueg.github.io/bixverse/reference/rs_gower_dist.md)
  **\[experimental\]** : Calculates the Gower distance for a given
  matrix
- [`rs_gse_geom_elim()`](https://gregorlueg.github.io/bixverse/reference/rs_gse_geom_elim.md)
  **\[experimental\]** : Run hypergeometric enrichment over the gene
  ontology
- [`rs_gse_geom_elim_list()`](https://gregorlueg.github.io/bixverse/reference/rs_gse_geom_elim_list.md)
  **\[experimental\]** : Run hypergeometric enrichment a list of target
  genes over the gene ontology
- [`rs_gsva()`](https://gregorlueg.github.io/bixverse/reference/rs_gsva.md)
  **\[experimental\]** : Rust version of the GSVA algorithm
- [`rs_h5ad_data()`](https://gregorlueg.github.io/bixverse/reference/rs_h5ad_data.md)
  **\[experimental\]** : Load in h5ad data via Rust
- [`rs_hamming_dist()`](https://gregorlueg.github.io/bixverse/reference/rs_hamming_dist.md)
  **\[experimental\]** : Calculates the Hamming distance between
  categorical columns
- [`rs_harmony()`](https://gregorlueg.github.io/bixverse/reference/rs_harmony.md)
  **\[experimental\]** : Harmony batch correction in Rust
- [`rs_harmony_v2()`](https://gregorlueg.github.io/bixverse/reference/rs_harmony_v2.md)
  **\[experimental\]** : Harmony batch correction in Rust (version 2)
- [`rs_hedges_g()`](https://gregorlueg.github.io/bixverse/reference/rs_hedges_g.md)
  **\[experimental\]** : Calculate the Hedge's G effect
- [`rs_hypergeom_test()`](https://gregorlueg.github.io/bixverse/reference/rs_hypergeom_test.md)
  **\[experimental\]** : Run a single hypergeometric test.
- [`rs_hypergeom_test_list()`](https://gregorlueg.github.io/bixverse/reference/rs_hypergeom_test_list.md)
  **\[experimental\]** : Run a hypergeometric test over a list of target
  genes
- [`rs_ica_iters()`](https://gregorlueg.github.io/bixverse/reference/rs_ica_iters.md)
  **\[experimental\]** : Run ICA over a given no_comp with random
  initilisations of w_init
- [`rs_ica_iters_cv()`](https://gregorlueg.github.io/bixverse/reference/rs_ica_iters_cv.md)
  **\[experimental\]** : Run ICA with cross-validation and random
  initialsiation
- [`rs_jaccard_row_integers()`](https://gregorlueg.github.io/bixverse/reference/rs_jaccard_row_integers.md)
  **\[experimental\]** : Calculate rapidbly Jaccard similarities between
  rows
- [`rs_kbet()`](https://gregorlueg.github.io/bixverse/reference/rs_kbet.md)
  **\[experimental\]** : Calculate kBET type scores
- [`rs_knn_label_propagation()`](https://gregorlueg.github.io/bixverse/reference/rs_knn_label_propagation.md)
  **\[experimental\]** : kNN label propagation
- [`rs_knn_mat_to_edge_list()`](https://gregorlueg.github.io/bixverse/reference/rs_knn_mat_to_edge_list.md)
  **\[experimental\]** : Flatten kNN matrix to edge list
- [`rs_generate_ligand_target_influence()`](https://gregorlueg.github.io/bixverse/reference/rs_generate_ligand_target_influence.md)
  **\[experimental\]** : Generate the ligand to target influence
  matrices
- [`rs_ligand_activity_scores()`](https://gregorlueg.github.io/bixverse/reference/rs_ligand_activity_scores.md)
  **\[experimental\]** : Calculate the NicheNet ligand activity scores
- [`rs_mad_outlier()`](https://gregorlueg.github.io/bixverse/reference/rs_mad_outlier.md)
  **\[experimental\]** : Calculate MAD outlier detection in Rust.
- [`rs_mc_aucell()`](https://gregorlueg.github.io/bixverse/reference/rs_mc_aucell.md)
  **\[experimental\]** : Calculate AUCell in Rust (for meta cells)
- [`rs_mc_hvg()`](https://gregorlueg.github.io/bixverse/reference/rs_mc_hvg.md)
  **\[experimental\]** : Meta cells highly variable genes
- [`rs_mc_pca()`](https://gregorlueg.github.io/bixverse/reference/rs_mc_pca.md)
  **\[experimental\]** : PCA on MetaCells (sparse data)
- [`rs_mc_scenic()`](https://gregorlueg.github.io/bixverse/reference/rs_mc_scenic.md)
  **\[experimental\]** : SCENIC on MetaCells
- [`rs_mitch_calc()`](https://gregorlueg.github.io/bixverse/reference/rs_mitch_calc.md)
  : Calculate mitch enrichment leveraging Rust under the hood
- [`rs_mutual_info()`](https://gregorlueg.github.io/bixverse/reference/rs_mutual_info.md)
  **\[experimental\]** : Calculates the mutual information matrix
- [`rs_nmf_multi_mc()`](https://gregorlueg.github.io/bixverse/reference/rs_nmf_multi_mc.md)
  **\[experimental\]** : Run multiple NMF (HALS) restarts on MetaCells
- [`rs_nmf_single_mc()`](https://gregorlueg.github.io/bixverse/reference/rs_nmf_single_mc.md)
  **\[experimental\]** : Run NMF (HALS) on MetaCells
- [`rs_onto_semantic_sim()`](https://gregorlueg.github.io/bixverse/reference/rs_onto_semantic_sim.md)
  **\[experimental\]** : Calculate the semantic similarity in an
  ontology
- [`rs_onto_semantic_sim_mat()`](https://gregorlueg.github.io/bixverse/reference/rs_onto_semantic_sim_mat.md)
  **\[experimental\]** : Calculate the semantic similarity in an
  ontology
- [`rs_onto_sim_wang()`](https://gregorlueg.github.io/bixverse/reference/rs_onto_sim_wang.md)
  **\[experimental\]** : Calculate the Wang similarity for specific
  terms
- [`rs_onto_sim_wang_mat()`](https://gregorlueg.github.io/bixverse/reference/rs_onto_sim_wang_mat.md)
  **\[experimental\]** : Calculate the Wang similarity matrix for an
  ontology
- [`rs_ot_harmonic_sum()`](https://gregorlueg.github.io/bixverse/reference/rs_ot_harmonic_sum.md)
  **\[experimental\]** : Calculate the OT harmonic sum
- [`rs_page_rank()`](https://gregorlueg.github.io/bixverse/reference/rs_page_rank.md)
  **\[experimental\]** : Rust version of calcaluting the personalised
  page rank
- [`rs_page_rank_parallel()`](https://gregorlueg.github.io/bixverse/reference/rs_page_rank_parallel.md)
  **\[experimental\]** : Calculate massively parallelised personalised
  page rank scores
- [`rs_phyper()`](https://gregorlueg.github.io/bixverse/reference/rs_phyper.md)
  **\[experimental\]** : Calculate the hypergeometric test in Rust
- [`rs_pointwise_mutual_info()`](https://gregorlueg.github.io/bixverse/reference/rs_pointwise_mutual_info.md)
  **\[experimental\]** : Calculates the point wise mutual information
- [`rs_prcomp()`](https://gregorlueg.github.io/bixverse/reference/rs_prcomp.md)
  **\[experimental\]** : Rust implementation of prcomp
- [`rs_prepare_whitening()`](https://gregorlueg.github.io/bixverse/reference/rs_prepare_whitening.md)
  **\[experimental\]** : Prepare the data for whitening
- [`rs_random_svd()`](https://gregorlueg.github.io/bixverse/reference/rs_random_svd.md)
  **\[experimental\]** : Run randomised SVD over a matrix
- [`rs_range_norm()`](https://gregorlueg.github.io/bixverse/reference/rs_range_norm.md)
  **\[experimental\]** : Apply a range normalisation on a vector.
- [`rs_rank_matrix_col()`](https://gregorlueg.github.io/bixverse/reference/rs_rank_matrix_col.md)
  **\[experimental\]** : Gene rank matrix
- [`rs_rank_matrix_col_stable()`](https://gregorlueg.github.io/bixverse/reference/rs_rank_matrix_col_stable.md)
  **\[experimental\]** : Stable-gene rank matrix
- [`rs_rbf_function()`](https://gregorlueg.github.io/bixverse/reference/rs_rbf_function.md)
  **\[experimental\]** : Apply a Radial Basis Function
- [`rs_rbf_function_mat()`](https://gregorlueg.github.io/bixverse/reference/rs_rbf_function_mat.md)
  **\[experimental\]** : Apply a Radial Basis Function (to a matrix)
- [`rs_rbh_cor()`](https://gregorlueg.github.io/bixverse/reference/rs_rbh_cor.md)
  **\[experimental\]** : Generate reciprocal best hits based on
  correlations
- [`rs_rbh_sets()`](https://gregorlueg.github.io/bixverse/reference/rs_rbh_sets.md)
  **\[experimental\]** : Generate reciprocal best hits based on set
  similarities
- [`rs_sc_knn()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_knn.md)
  **\[experimental\]** : Generates the kNN graph
- [`rs_sc_knn_w_dist()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_knn_w_dist.md)
  **\[experimental\]** : Generates the kNN graph with additional
  distances
- [`rs_sc_snn()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_snn.md)
  **\[experimental\]** : Generates the sNN graph for igraph
- [`rs_set_similarity()`](https://gregorlueg.github.io/bixverse/reference/rs_set_similarity.md)
  **\[experimental\]** : Set similarities
- [`rs_set_similarity_list()`](https://gregorlueg.github.io/bixverse/reference/rs_set_similarity_list.md)
  **\[experimental\]** : Set similarities over one list
- [`rs_set_similarity_list2()`](https://gregorlueg.github.io/bixverse/reference/rs_set_similarity_list2.md)
  **\[experimental\]** : Set similarities over two list
- [`rs_singscore_multi()`](https://gregorlueg.github.io/bixverse/reference/rs_singscore_multi.md)
  **\[experimental\]** : Rust version of singscore for many gene sets
- [`rs_singscore_permutation_test()`](https://gregorlueg.github.io/bixverse/reference/rs_singscore_permutation_test.md)
  **\[experimental\]** : Rust version of the singscore permutation test
- [`rs_singscore_single()`](https://gregorlueg.github.io/bixverse/reference/rs_singscore_single.md)
  **\[experimental\]** : Rust version of singscore for a single gene set
- [`rs_snf()`](https://gregorlueg.github.io/bixverse/reference/rs_snf.md)
  **\[experimental\]** : Similarity network fusion
- [`rs_sparse_dict_dgrdl()`](https://gregorlueg.github.io/bixverse/reference/rs_sparse_dict_dgrdl.md)
  **\[experimental\]** : Generate a sparse dictionary with DGRDL
- [`rs_sparse_dict_dgrdl_grid_search()`](https://gregorlueg.github.io/bixverse/reference/rs_sparse_dict_dgrdl_grid_search.md)
  **\[experimental\]** : Generate a sparse dictionary with DGRDL
- [`rs_spectral_clustering()`](https://gregorlueg.github.io/bixverse/reference/rs_spectral_clustering.md)
  **\[experimental\]** : Rust implementation of spectral clustering
- [`rs_spectral_clustering_sim()`](https://gregorlueg.github.io/bixverse/reference/rs_spectral_clustering_sim.md)
  **\[experimental\]** : Rust implementation of spectral clustering
- [`rs_ssgsea()`](https://gregorlueg.github.io/bixverse/reference/rs_ssgsea.md)
  **\[experimental\]** : Rust version of the ssGSEA algorithm
- [`rs_tied_diffusion_parallel()`](https://gregorlueg.github.io/bixverse/reference/rs_tied_diffusion_parallel.md)
  **\[experimental\]** : Calculate massively parallelised tied diffusion
  scores
- [`rs_tom()`](https://gregorlueg.github.io/bixverse/reference/rs_tom.md)
  **\[experimental\]** : Calculates the TOM over an affinity matrix
- [`rs_upper_triangle_to_dense()`](https://gregorlueg.github.io/bixverse/reference/rs_upper_triangle_to_dense.md)
  **\[experimental\]** : Reconstruct a matrix from a flattened upper
  triangle vector
- [`rs_upper_triangle_to_sparse()`](https://gregorlueg.github.io/bixverse/reference/rs_upper_triangle_to_sparse.md)
  **\[experimental\]** : Generate sparse data from an upper triangle
- [`rs_wnn()`](https://gregorlueg.github.io/bixverse/reference/rs_wnn.md)
  **\[experimental\]** : Run the weighted nearest neighbour algorithm
