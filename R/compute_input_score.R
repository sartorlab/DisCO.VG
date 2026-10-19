#' @title Compute Input Score Module
#'
#' @description
#' DisCO-VG input score functionality. Calculates for subscores and
#' integrates these subscore together: GWAS, distance, eQTL, EnTDef.
#'
#' @importFrom foreach foreach `%dopar%`
#' @importFrom doParallel registerDoParallel
#' @importFrom GenomicRanges distanceToNearest findOverlaps makeGRangesFromDataFrame
#' @importFrom Matrix sparseMatrix summary
#' @importFrom pracma sigmoid
#' @import chipenrich.data
#'
#' @name compute_input_score
NULL

#' Compute DisCO-VG input scores for each variant-gene pair
#'
#' Compute DisCO-VG input scores for each relevant variant-gene pair.
#' This function used GWAS variants and a gene list to integrate multiple lines
#' of functional evidence for either variant, gene, or variant-gene pair.
#' These scores are represented in a main matrix, "MVG".
#'
#' @param variants List of GWAS variants being analyzed. Typically a product of
#'    the `filter_variants` module.
#' @param gene_list List of genes being considered in the DisCO-VG study.
#'    Defaults to built-in gene list that is filtered for genes that have at least one GO
#'    BP (Biological Process) term.
#' @param df_eQTL eQTL summary statistics for calculation of the eQTL score.
#'    Should be prefiltered to user specifications (e.g., p < 0.05).
#' @param subscore_weights Weights to assign each main matrix subscore; default is 1:3:1:1
#'    for GWAS, distance, eQTL, and EnTDef scores.
#' @param genome Genome reference. Only supports hg19 and hg38 right now.
#' @param n_cores Number of cores to allocate for computing distance scores.
#'
#' @return sparse matrix containing integrated, normalized input scores for each variant-pair.
#'
#' @export
compute_input_score <- function(variants,
                                gene_list = NULL,
                                df_eQTL,
                                subscore_weights = list(
                                  GWAS = 1,
                                  dist = 3,
                                  eQTL = 1,
                                  EnTDef = 1
                                ),
                                genome = c("hg19", "hg38"),
                                n_cores = 4) {
  # Load default gene_list if necessary
  if (is.null(gene_list)) {
    gene_list = load_gene_list()
  }

  # Verify that the data provided contains the expected columns
  stopifnot(
    is.data.frame(variants),
    all(c("SNP", "CHR", "BP", "P") %in% colnames(variants)),
    is.data.frame(gene_list),
    nrow(gene_list) > 2,
    all(c("ensembl_id", "entrez_id", "symbol") %in% colnames(gene_list)),
    is.data.frame(df_eQTL),
    all(c("SNP", "P") %in% colnames(df_eQTL)),
    any(c("ensembl_id", "entrez_id", "symbol") %in% colnames(df_eQTL))
  )

  # Compute distances between variant-gene pairs
  message("Computing variant-gene distances...")
  distance_in_bp = distance_matrix(variants, gene_list, genome, n_cores)

  # Construct list of matrices with subscores
  message("Computing subscores...")
  subscore_matrices = list()
  subscore_matrices[["GWAS"]] = compute_gwas_score(variants, gene_list)
  subscore_matrices[["dist"]] = compute_distance_score(distance_in_bp)
  subscore_matrices[["eQTL"]] = compute_eqtl_score(variants, gene_list, df_eQTL)
  subscore_matrices[["EnTDef"]] = compute_entdef_score(variants, gene_list, distance_in_bp, genome)

  # Integrate subscores into final score matrix
  message("Integrating subscores into final score matrix...")
  message(paste(intersect(names(subscore_matrices), names(subscore_weights)), collapse = "; "))
  input_mat = integrate_subscores(subscore_matrices, subscore_weights)

  return(input_mat)
}


#' Compute GWAS score for each variant
#'
#' Computes a subscore for each relevant variant from GWAS summary
#' statistics. Using the GWAS p-value, we transform that into a score between
#' 0-1 using a sigmoid function. Every gene linked to this given variant receives
#' the same GWAS subscore.
#'
#' @param variants List of GWAS variants being analyzed. Typically a product of
#'    the `filter_variants` module. Expected columns:
#'     - SNP
#'     - GWAS significance label (from the `label_variants` function)
#' @param gene_list List of genes being considered in the DisCO-VG study.
#'    Defaults to built-in gene list that is filtered for genes that have at least one GO
#'    BP (Biological Process) term.
#' @param sigmoid_params Parameters for the sigmoid function, implemented by the
#'    `pracma::sigmoid` function.
#' @return Sparse matrix with GWAS subscores for all relevant variant-gene pairs.
#' @export
compute_gwas_score <- function(variants,
                               gene_list,
                               sigmoid_params = list(
                                 a = 0.5,
                                 b = 10
                               )) {
  message("Computing GWAS score...")
  gwas_sigmoid_scores = pracma::sigmoid(-log10(variants$P),
                                        sigmoid_params$a,
                                        sigmoid_params$b)

  idx_nonsig = grepl("nonsignificant", variants$label)
  gwas_sigmoid_scores[idx_nonsig] = 0 # set (proximal) nonsignificant scores to 0

  nr = nrow(variants) # number of rows
  nc = nrow(gene_list) # number of columns
  nonzero_indices = which(gwas_sigmoid_scores > 0)
  x = gwas_sigmoid_scores[nonzero_indices]

  mat_GWAS = Matrix::sparseMatrix(
    i = rep.int(nonzero_indices, times = nc),
    j = rep(seq_len(nc), each = length(nonzero_indices)),
    x = rep(x, times = nc),
    dims = c(nr, nc),
    dimnames = list(variants$SNP, gene_list$ensembl_id)
  )
  return(mat_GWAS)
}

#' Compute distance score for each relevant variant-gene pair
#'
#' Whenever applicable, compute our heuristic distance score as part of the
#' DisCO-VG input score. Any variant within 1Mb of a given gene gets a score
#' between 0 and 1.
#'
#' In the current version, the scores are:
#'  * 1.0 if < 1kb
#'  * 0.8 if between 1-5kb
#'  * 0.6 if between 5-20kb
#'  * 0.4 if between 20-100kb
#'  * 0.2 if between 100-300kb
#'  * 0.1 if between 300kb-1Mb
#'  * 0 o.w.
#'
#' @param mat_dist.bp Matrix of variant-gene distances in base pairs, as computed by the `distance_matrix()` function
#' @return Sparse matrix of distance scores for each variant-gene pair
#' @export
compute_distance_score <- function(mat_dist.bp) {
  message("Computing distance score...")
  df_dist = Matrix::summary(mat_dist.bp)

  df_dist = dplyr::mutate(df_dist,
                          score = dplyr::case_when(
                            x < 1e3 ~ 1.0,
                            x < 5e3 ~ 0.8,
                            x < 20e3 ~ 0.6,
                            x < 100e3 ~ 0.4,
                            x < 300e3 ~ 0.2,
                            x < 1e6 ~ 0.1,
                            .default = 0
                          ))

  df_dist = dplyr::filter(df_dist, score > 0)

  mat_dist.score = Matrix::sparseMatrix(
    i = df_dist$i,
    j = df_dist$j,
    x = df_dist$score,
    dims = c(nrow(mat_dist.bp), ncol(mat_dist.bp)),
    dimnames = mat_dist.bp@Dimnames
  )

  return(mat_dist.score)
}

#' Compute eQTL score for each variant-gene pair
#'
#' Computes a subscore for each relevant variant-gene pair from eQTL summary
#' statistics. Using the eQTL p-value, we transform that into a score between
#' 0-1 using a sigmoid function. If there are duplicate SNP-gene pairs, only
#' the most significant eQTL is kept.
#'
#' @param variants List of GWAS variants being analyzed. Typically a product of
#'    the `filter_variants` module.
#' @param gene_list List of genes being considered in the DisCO-VG study.
#'    Defaults to built-in gene list that is filtered for genes that have at least one GO
#'    BP (Biological Process) term.
#' @param df_eQTL eQTL summary statistics for calculation of the eQTL score.
#'    Should be prefiltered to user specifications (e.g., p < 0.05).
#' @param eQTL_sigmoid_params Parameters for the sigmoid function, implemented
#'    by the `pracma::sigmoid` function.
#' @return Sparse matrix with eQTL subscores for all relevant variant-gene pairs.
#' @export
compute_eqtl_score <- function(variants,
                               gene_list,
                               df_eQTL,
                               eQTL_sigmoid_params = list(
                                 a = 0.5,
                                 b = 2
                               )) {
  message("Computing eQTL score...")
  df_eQTL.filtered = dplyr::filter(df_eQTL,
                                   SNP %in% variants$SNP)
  df_eQTL.filtered$eQTL_score = pracma::sigmoid(-log10(df_eQTL.filtered$P),
                                                a = eQTL_sigmoid_params$a,
                                                b = eQTL_sigmoid_params$b)

  # Map to EnsemblID if necessary
  if ("ensembl_id" %in% colnames(df_eQTL.filtered)) {
    # First priority, but skip to post if statement
  } else if ("symbol" %in% colnames(df_eQTL.filtered)) {
    df_eQTL.filtered = dplyr::inner_join(df_eQTL.filtered,
                                         gene_list[, c("symbol", "ensembl_id")],
                                         by = "symbol",
                                         relationship = "many-to-many")
  } else {
    # assuming ENTREZID is the only remaining option
    df_eQTL.filtered = dplyr::inner_join(df_eQTL.filtered,
                                         gene_list[, c("entrez_id", "ensembl_id")],
                                         by = "entrez_id",
                                         relationship = "many-to-many")
  }

  # If applicable, filter out SNP, ENSEMBL pair duplicates
  df_eQTL.filtered <- df_eQTL.filtered %>%
    group_by(SNP, ensembl_id) %>%
    slice_max(eQTL_score) %>%
    ungroup()

  mat_eQTL = Matrix::sparseMatrix(
      i = match(df_eQTL.filtered$SNP, variants$SNP),
      j = match(df_eQTL.filtered$ensembl_id, gene_list$ensembl_id),
      x = df_eQTL.filtered$eQTL_score,
      dims = c(nrow(variants), nrow(gene_list)),
      dimnames = list(variants$SNP, gene_list$ensembl_id)
  )
  return(mat_eQTL)
}

#' Compute EnTDef subscore for DisCO-VG
#'
#' Computes EnTDef (Enhancer-Target gene Definition) score for a given
#' variant-gene pair. For any variant in our cell-type-agnostic enhancer-target gene list,
#' we assign a score of 1 to the gene in the EnTDef list and 0.6 if it is not
#' the matching gene but is within 1Mb of the variant.
#'
#' @param variants List of variants being included in the DisCO-VG analysis,
#'    with at minimum the column `SNP`, `CHR`, `BP`
#' @param gene_list List of genes being considered in the DisCO-VG study.
#'  Defaults to built-in gene list that is filtered for genes that have at least one GO
#'  BP (Biological Process) term.
#' @param mat_dist.bp Matrix of variant-gene distances in base pairs as calculated by the `distance_matrix()` function.
#' @param genome Reference genome for EnTDef list
#' @return Sparse matrix of variant-gene pairs where the entries reflect EnTDef scores.
#' @export
compute_entdef_score <- function(variants,
                                 gene_list,
                                 mat_dist.bp,
                                 genome = c("hg19", "hg38")) {
  message("Computing EnTDef score...")
  genome = match.arg(genome)
  entdef = load_entdef(genome)

  # Define scores here
  EnTDef_mismatch_score = 0.6
  EnTDef_max_score = 1.0

  df_dist = Matrix::summary(mat_dist.bp)

  # Identify which variants are in an enhancer
  message("Filtering for variants in enhancers...")
  gr_var = GenomicRanges::makeGRangesFromDataFrame(variants,
                                                   start.field = "BP",
                                                   end.field = "BP")
  gr_entdef = GenomicRanges::makeGRangesFromDataFrame(entdef, keep.extra.columns = T)
  overlaps = GenomicRanges::findOverlaps(gr_var, gr_entdef)

  entdef.filtered = entdef[overlaps@to,]
  entdef.filtered$var_in_enh = variants$SNP[overlaps@from]
  entdef.filtered$var_in_enh.idx = match(entdef.filtered$var_in_enh, rownames(mat_dist.bp))
  entdef.filtered = dplyr::left_join(entdef.filtered, gene_list,
                                     by = dplyr::join_by(gene_id == entrez_id),
                                     relationship = "many-to-many")
  # account for NAs, even though this shouldn't be possible by design in variants
  # (but expected in genes)
  entdef.filtered = dplyr::filter(entdef.filtered,
                                  !is.na(ensembl_id), !is.na(var_in_enh.idx))
  entdef.filtered$gene_idx = match(entdef.filtered$ensembl_id, colnames(mat_dist.bp))

  # Assign every gene within 1Mb of these variants the "mismatch score"
  message("Assigning EnTDef scores to VEG triplets within 1Mb of each other...")
  df_entdef = dplyr::filter(df_dist,
                            i %in% entdef.filtered$var_in_enh.idx,
                            x < 1e6)

  # Assign every variant-gene pair found in the EnTDef list the "max score"
  message("Assigning exact matches of VG pairs in the EnTDef list...")
  # Filter for NA genes, just to cover edge cases
  entdef.filtered = dplyr::filter(entdef.filtered,
                                  !is.na(gene_idx))
  df_entdef = dplyr::left_join(df_entdef,
                               entdef.filtered,
                               by = dplyr::join_by(i == var_in_enh.idx,
                                                   j == gene_idx))

  df_entdef = dplyr::mutate(df_entdef,
                            EnTDef = dplyr::if_else(is.na(gene_id),
                                                    EnTDef_mismatch_score,
                                                    EnTDef_max_score))

  # Set up sparse matrix
  mat_EnTDef = Matrix::sparseMatrix(
    i = df_entdef$i,
    j = df_entdef$j,
    x = df_entdef$EnTDef,
    dims = dim(mat_dist.bp),
    dimnames = mat_dist.bp@Dimnames
  )

  return(mat_EnTDef)
}

#' Integrate subscores into main matrix input scores
#'
#' Given a list of matrices containing subscores and a list with their respective
#' weights, add the scores together and min-max normalize the scores. Additionally,
#' if distance and/or EnTDef scores are a part of the input, only variant-gene
#' pairs within 1Mb of each other and/or in an EnTDef are considered in the final
#' input matrix.
#'
#' This function has been generalized to support linear combinations of any kinds
#' of subscores, not just the ones that this package offers. Only subscores for
#' which both the (sparse) matrix and the weight are provided will be integrated
#' into the final input score matrix.
#'
#' @param matrices List of named matrices containing subscores for the same
#'    variant-gene pairs in each entry, i.e., all matrices' dimensions and order
#'    of variants and genes in rows/columns respectively should be the same.
#' @param subscore_weights List of named weights for the computation of the linear
#'    combination of the subscore matrices.
#' @return MVG: Sparse matrix with DisCO-VG input scores for each relevant variant-gene
#'    pair
#' @export
integrate_subscores <- function(matrices,
                                subscore_weights) {
  # verify the matrices are all the same dimensions
  stopifnot(all(sapply(matrices[-1], dim) == dim(matrices[[1]])))

  # verify the weights and matrices have the same keys
  stopifnot(all(names(matrices) %in% names(subscore_weights)))

  # iteratively add weighted matrices together
  shared_keys = intersect(names(matrices), names(subscore_weights))
  stopifnot(length(shared_keys) > 0)

  integrated_mat = matrices[[shared_keys[1]]] * subscore_weights[[shared_keys[1]]]
  if (length(shared_keys) > 1) {
    for (i in 2:length(shared_keys)) {
      integrated_mat = integrated_mat + matrices[[shared_keys[i]]] * subscore_weights[[shared_keys[i]]]
    }
  }

  # Limit input scores to variant-gene pairs that have either a nonzero distance score
  # (i.e., the variant is within 1Mb of the gene) or a nonzero EnTDef score
  # (i.e., the variant overlaps an enhancer linked to the gene)
  if (all(c("dist", "EnTDef") %in% shared_keys)) {
    integrated_mat = integrated_mat * (matrices[["dist"]] + matrices[["EnTDef"]] > 0)
  } else if ("dist" %in% shared_keys) {
    integrated_mat = integrated_mat * (matrices[["dist"]] > 0)
  } else if ("EnTDef" %in% shared_keys) {
    integrated_mat = integrated_mat * (matrices[["EnTDef"]] > 0)
  }

  return(normalize(integrated_mat))
}

#' Generate pool of nonzero variant-gene scores to sample from for bagging
#'
#' This function takes in a DisCO-VG input score matrix (MVG, where the rows
#' represent variants and the columns represent genes) and extracts the nonzero
#' entries from this matrix. It stratifies these variant-gene pairs with nonzero
#' scores into three groups:
#'  - `positive`: GWAS-significant and -moderate variants
#'  - `prox_nonsig`: one of our negative controls; variants labeled "proximal
#'    nonsignificant" in your data
#'  - `nonsig`: second set of negative controls; variants labeled "nonsignificant"
#'    in your data
#'
#' @param variants List of GWAS variants being analyzed. Typically a product of
#'    the `filter_variants` module.
#' @param MVG Sparse matrix with DisCO-VG input scores for each relevant
#'    variant-gene pair.
#' @return List of three dataframes containing variant-gene pairs with non-zero scores:
#'  - `positive`
#'  - `prox_nonsig`
#'  - `nonsig`
#' @export
generate_sample_pool <- function(variants, MVG) {
  sample_pool_VG = list()

  nonzero_entries = Matrix::summary(MVG)
  nonzero_entries$variant = rownames(MVG)[nonzero_entries$i]
  nonzero_entries$gene = colnames(MVG)[nonzero_entries$j]
  nonzero_entries = dplyr::rename(nonzero_entries, input_score = x)

  positive_vars = variants$SNP[variants$label %in% c("significant", "moderate")]
  prox_nonsig_vars = variants$SNP[variants$label == "proximal nonsignificant"]
  nonsig_vars = variants$SNP[variants$label == "nonsignificant"]

  sample_pool_VG[["positive"]] = dplyr::filter(nonzero_entries, variant %in% positive_vars)
  sample_pool_VG[["prox_nonsig"]] = dplyr::filter(nonzero_entries, variant %in% prox_nonsig_vars)
  sample_pool_VG[["nonsig"]] = dplyr::filter(nonzero_entries, variant %in% nonsig_vars)

  return(sample_pool_VG)
}


#' Compute variant-gene distance matrix
#'
#' Compute genomic distance between variant and gene for each variant-gene pair
#' being studied. This function serves to make downstream tasks easier, but
#' it can also be used by users designing their own input scores.
#'
#' @param variants List of variants with at minimum columns `SNP`, `CHR`, `BP`.
#' @param gene_list List of genes in study. Should be in a similar format as
#'   `default_gene_list`.
#' @param genome Reference genome; currently supports hg19 and hg38.
#' @param n_cores Number of cores; recommended to use multiple to compute
#'  distance matrix more efficiently.
#' @return Sparse matrix with variants as rows and genes as columns, where the entries
#'   are the variant-gene distance (in bp), if they are on the same chromosome,
#'   and 0 otherwise. If the distance between variant gene is zero, it gets converted into a pseudocount of 1e-16.
#' @export
distance_matrix <- function(variants,
                            gene_list = NULL,
                            genome = c("hg19", "hg38"),
                            n_cores) {
  # Following Raymond's example, check for valid genome arguments
  genome = match.arg(genome)

  if (is.null(gene_list)) {
    gene_list = load_gene_list()
  }

  gr_vars = GenomicRanges::makeGRangesFromDataFrame(variants,
                                                    start.field = "BP",
                                                    end.field = "BP")
  # # Initialize dataframe that will be passed on to sparseMatrix
  # mat_values = data.frame(i = integer(),
  #                         j = integer(),
  #                         x = numeric())

  if (genome == "hg19") {
    gr_tss = chipenrich.data::tss.hg19
  } else if (genome == "hg38") {
    gr_tss = chipenrich.data::tss.hg38
  }

  registerDoParallel(n_cores)

  genes_in_chipenrich = which(def_gene_list$entrez_id %in% gr_tss$gene_id)

  mat_list = foreach(entrez_id = gene_list$entrez_id[genes_in_chipenrich]) %dopar% {
    vg_distances = GenomicRanges::distanceToNearest(gr_vars, gr_tss[gr_tss$gene_id == entrez_id])

    gene_indices = which(gene_list$entrez_id == entrez_id)
    nr = length(vg_distances@from)
    nc = length(gene_indices)
    data.frame(i = rep(vg_distances@from, times = nc),
               j = rep(gene_indices, each = nr),
               x = vg_distances@elementMetadata@listData[["distance"]])
  }

  mat_values = data.table::rbindlist(mat_list)

  # pseudo count any VG distances that are zero to store data in sparse matrix
  mat_values$x[mat_values$x < 1] = 1e-16

  mat_dist.bp = Matrix::sparseMatrix(
    i = mat_values$i,
    j = mat_values$j,
    x = mat_values$x,
    dims = c(nrow(variants), nrow(gene_list)),
    dimnames = list(variants$SNP, gene_list$ensembl_id)
  )

  return(mat_dist.bp)
}

# Helper function to min-max normalize
normalize <- function(x) {
  min_x = min(x,na.rm = TRUE)
  max_x = max(x,na.rm = TRUE)
  return((x - min_x) / (max_x-min_x))
}

#' Load default gene list
#'
#' Loads default gene list if no gene list (dataframe of genes with Entrez ID,
#' EnsemblID, and HGNC gene symbol) is supplied by the user. This default list
#' contains each gene that has at least one GO Biological Process (BP) term
#' for computing gene-gene similarity scores.
#'
#' @return Default gene list
load_gene_list <- function() {
  e <- new.env(parent = emptyenv())
  load(system.file("data", "default_gene_list.rda", package = "DisCO.VG"),
       envir = e)
  return(e$default_gene_list)
}

#' Load EnTDef list
#'
#' Loads default gene list if no gene list (dataframe of genes with Entrez ID,
#' EnsemblID, and HGNC gene symbol) is supplied by the user. This default list
#' contains each gene that has at least one GO Biological Process (BP) term
#' for computing gene-gene similarity scores.
#'
#' @return Default gene list
load_entdef <- function(genome = c("hg19", "hg38")) {
  genome = match.arg(genome)

  e <- new.env(parent = emptyenv())
  if (genome == "hg19") {
    load(system.file("data", "entdef.hg19.rda", package = "DisCO.VG"),
         envir = e)
    return(e$entdef.hg19)
  } else if (genome == "hg38") {
    load(system.file("data", "entdef.hg38.rda", package = "DisCO.VG"),
         envir = e)
    return(e$entdef.hg38)
  }
}
