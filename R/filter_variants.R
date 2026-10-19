#' @title Filter for GWAS variants with functional evidence
#'
#' @description
#' Code to filter GWAS variants for those in ATAC-Seq peaks.
#' Also labels variants based on GWAS significance and prunes superfluous nonsignificant variants.
#'
#' @importFrom data.table fread
#' @importFrom dplyr case_when if_else filter mutate pull
#' @importFrom GenomicRanges findOverlaps makeGRangesFromDataFrame
#' @importFrom magrittr `%>%`
#'
#' @name filter_variants
NULL

# --- 1. Filter variants ---
#' Filter for GWAS variants in ATAC-Seq peaks
#'
#' Filter for GWAS variants that are in ATAC-Seq peaks. This function uses standard GRanges functionality.
#'
#' @param path_GWAS Path to GWAS summary statistics file, with at minimum columns "SNP", "CHR", "BP", and "P".
#' @param path_ATAC_regions Path to ATAC-Seq peaks file, with at minimum columns "chr", "start", "end", and "peak_ID".
#' @return A dataframe of GWAS variants annotated with ATAC-Seq peak IDs
#' @export
filter_variants_in_ATAC <- function(path_GWAS, path_ATAC_regions) {
  # Load data
  message("Loading GWAS data...")
  df_gwas = data.table::fread(path_GWAS, data.table = FALSE)
  message("Loading ATAC-Seq data...")
  df_atac = data.table::fread(path_ATAC_regions, data.table = FALSE)

  # check if files are empty; if so, throw error
  if (nrow(df_gwas) == 0 | nrow(df_atac) == 0) {
    stop("One or both of the input files are empty. Please check your input files.")
  }

  # Convert p-value to numeric if necessary
  if (!is.numeric(df_gwas$P[1])) {
    dplyr::mutate(df_gwas,
                  P = as.numeric(P))
  }

  # validate data (columns are what we expect)
  # NOTE: as of this version, we are being strict with column names
  # GWAS
  if (!all(c("SNP", "CHR", "BP", "P") %in% colnames(df_gwas))) {
    stop("Expected Munge-standardized GWAS column names 'SNP', 'CHR', 'BP', and 'P'.
    Provided GWAS data is missing necessary columns.
    Please ensure the GWAS summary statistics file has these specified columns.")
  }

  # ATAC-Seq
  if (!("peak_ID" %in% colnames(df_atac))) {
    stop("Provided ATAC-Seq data is missing necessary column 'peak_ID'.
    Please ensure the ATAC-Seq peaks file has a column named 'peak_ID' that contains unique identifiers for each peak.")
  }

  # TODO: check if variant ID column exists; otherwise, create one
  # NOTE: will implement later when I'm less strict about the column names.

  # Make GRanges objects for both datasets
  message("Filtering for GWAS variants in ATAC-Seq peaks...")
  gr_gwas = GenomicRanges::makeGRangesFromDataFrame(df_gwas, seqnames.field = "CHR", start.field = "BP", end.field = "BP")

  tryCatch(expr = {gr_atac = GenomicRanges::makeGRangesFromDataFrame(df_atac, keep.extra.columns = TRUE)},
           error = function(e) {
             msg = "Error in creating GRanges object for ATAC-Seq data. Please check your ATAC-Seq peak locations file.\nEnsure that the file has columns that the GenomicRanges package expects: 'chr', 'start', and 'end'."
             e$message = msg
             stop(e)
  })

  # Find overlaps between GWAS variants and ATAC-Seq peaks
  overlaps = GenomicRanges::findOverlaps(gr_gwas, gr_atac)

  # Keep only GWAS variants in ATAC-Seq peaks; annotate with ATAC peak ID
  df_gwas.filtered = df_gwas[overlaps@from,] %>%
    dplyr::mutate(peak_ID = gr_atac$peak_ID[overlaps@to])

  # if df_gwas.filtered is empty at this stage, throw warning (no need to throw error at this stage).
  if (nrow(df_gwas.filtered) == 0) {
    warning("No GWAS variants overlap with ATAC-Seq peaks.
    Please check your input files.
    Returning an empty dataframe.")
  }

  return(df_gwas.filtered)
}

# --- 2. Label variants ---
#' Label variants based on GWAS Significance
#'
#' This function uses the filtered variants from the previous step and applies labeling criteria
#' to categorize variants into "significant", "moderate", "proximal nonsignificant", and "nonsignificant" based on their p-values and proximity to significant variants.
#'
#' NOTE: This function was designed to run on a filtered set of GWAS variants
#' and thus may not scale to a full GWAS summary statistics file.
#'
#' @param df_gwas.filtered Dataframe of filtered GWAS variants annotated with ATAC-Seq peak IDs, ` function.output from the `filter_variants_in_ATAC` function.
#' @param dist_prox_nonsig Distance threshold for labeling a variant as "proximal nonsignificant". Default is 50kb.
#' @param GWAS.sig_threshold Threshold for labeling a variant as "significant". Default is 5e-8.
#' @param GWAS.moderate_threshold Threshold for labeling a variant as "moderate". Default is 1e-5.
#' @param GWAS.nonsig_threshold Threshold for labeling a variant as "nonsignificant". Default is 0.5.
#' @return A dataframe of GWAS variants annotated with GWAS significance labels ("significant", "moderate", "proximal nonsignificant", "nonsignificant") and ATAC-Seq peak IDs
#' @export
label_variants <- function(df_gwas.filtered, dist_prox_nonsig = 50000, GWAS.sig_threshold = 5e-8, GWAS.moderate_threshold = 1e-5, GWAS.nonsig_threshold = 0.5) {
  # Check if data is empty
  if (nrow(df_gwas.filtered) == 0) {
    stop("`label_variants` received an empty dataframe.")
  }

  stopifnot(all(c("SNP", "CHR", "BP", "P") %in% colnames(df_gwas.filtered)))

  # Ensuring numeric values where expected
  if (!is.numeric(df_gwas.filtered$BP[1])) {
    df_gwas.filtered$BP = as.numeric(df_gwas.filtered$BP)
  }

  if (!is.numeric(df_gwas.filtered$P[1])) {
    df_gwas.filtered$P = as.numeric(df_gwas.filtered$P)
  }

  # Label variants based on GWAS significance
  df_gwas.filtered = df_gwas.filtered %>%
    dplyr::filter(P < GWAS.moderate_threshold | P > GWAS.nonsig_threshold) %>%
    dplyr::mutate(label = dplyr::case_when(
      P < GWAS.sig_threshold ~ "significant",
      P < GWAS.moderate_threshold ~ "moderate",
      P > GWAS.nonsig_threshold ~ "nonsignificant"
    ))

  # For nonsignificant variants, check if they are within dist_prox_nonsig of a significant variant; if so, label as "proximal nonsignificant"
  significant_variants = df_gwas.filtered %>%
    dplyr::filter(label == "significant") %>%
    dplyr::mutate(start = BP - dist_prox_nonsig, end = BP + dist_prox_nonsig)
  nonsig_variants = df_gwas.filtered %>%
    dplyr::filter(label == "nonsignificant")

  # Use findOverlaps for a quicker way to label nonsignificant variants
  gr_sig = GenomicRanges::makeGRangesFromDataFrame(significant_variants)
  gr_nonsig = GenomicRanges::makeGRangesFromDataFrame(nonsig_variants,
                                                      start.field = "BP",
                                                      end.field = "BP")

  prox_nonsig_overlaps = GenomicRanges::findOverlaps(gr_nonsig, gr_sig)

  prox_nonsig_vars = nonsig_variants$SNP[unique(prox_nonsig_overlaps@from)]

  df_gwas.filtered = dplyr::mutate(df_gwas.filtered,
                                   label = dplyr::if_else(SNP %in% prox_nonsig_vars,
                                                          "proximal nonsignificant",
                                                          label))
  return(df_gwas.filtered)
}

# --- 3. Prune variants ---
#' Limit variant set for downstream analysis
#'
#' Prunes superfluous nonsignificant variants by labeling them as not included in the study (i.e. not in the sampling pool for bagging).
#' This function takes in the dataframe of labeled variants from the previous step and applies pruning criteria to determine which variants to keep
#' in the sampling pool for bagging.
#'
#' @param df_gwas.filtered.labeled A dataframe of GWAS variants that has
#'    1.) been filtered to only contain variants in ATAC-Seq peaks;
#'    2.) are labeled by GWAS significance.
#' @param n_variants Total number of variants being used during DisCO-VG model construction.
#'    Default is 15,000.
#' @param max_ratio_nonsig_to_sig Maximum ratio of nonsignificant variants to significant variants to include in the sampling pool.
#'    Default is 1.5 (i.e., max 1.5x more nonsignificant variants than significant variants).
#' @param seed Random seed for reproducibility. Default is 888.
#' @return A dataframe of GWAS variants annotated with GWAS significance labels
#'    ("significant", "moderate", "proximal nonsignificant", "nonsignificant"),
#'    ATAC-Seq peak IDs, and a boolean column indicating whether the variant is in the sampling pool for bagging.
#' @export
prune_variants <- function(df_gwas.filtered.labeled, n_variants = 15000, max_ratio_nonsig_to_sig = 1.5, seed = 888) {
  # Check if data is empty
  if (nrow(df_gwas.filtered.labeled) == 0) {
    stop("`prune_variants` received an empty dataframe.")
  }

  # Check if dataframe has necessary columns
  if (!all(c("SNP", "CHR", "BP", "P", "label") %in% colnames(df_gwas.filtered.labeled))) {
    stop("Dataframe is missing necessary columns. Please ensure the dataframe has columns 'SNP', 'CHR', 'BP', 'P', and 'label'.")
  }

  set.seed(seed)

  df_gwas.filtered.labeled$included_in_study = TRUE

  # Prune proximal nonsignificant variants if necessary
  n_sig = sum(df_gwas.filtered.labeled$label == "significant")
  n_prox_nonsig = sum(df_gwas.filtered.labeled$label == "proximal nonsignificant")
  if(n_prox_nonsig > max_ratio_nonsig_to_sig * n_sig) {
    n_to_sample = floor(max_ratio_nonsig_to_sig * n_sig)
    keep.prox_nonsig = sample(df_gwas.filtered.labeled %>% dplyr::filter(label == "proximal nonsignificant") %>% dplyr::pull(SNP), size = n_to_sample)
    df_gwas.filtered.labeled = df_gwas.filtered.labeled %>%
      dplyr::mutate(included_in_study = dplyr::if_else(label == "proximal nonsignificant" & !(SNP %in% keep.prox_nonsig),
                                                       FALSE,
                                                       included_in_study))
  }

  # Prune "redundant" nonsignficant variants by labeling as not included (in study);
  # all other variants are labeled TRUE.
  # TODO: edit this code to use "included_in_study" rather than just nrow and nonsig label
  if(sum(df_gwas.filtered.labeled$included_in_study) > n_variants) {
    # If we have more variants than n_variants, we need to prune some. We will keep all significant, moderate, and proximal nonsignificant variants,
    # and then sample from the rest to keep the total number at n_variants.
    n_to_sample = n_variants - sum((df_gwas.filtered.labeled$label != "nonsignificant") & (df_gwas.filtered.labeled$included_in_study))
    if (n_to_sample < 0) {
      warning("Number of significant, moderate, and proximal nonsignificant variants exceeds `n_variants``.
      Please consider increasing `n_variants` or adjusting GWAS thresholds as this will affect downstream analyses.
      This function will return more than `n_variants` in the sampling pool for bagging.
      Please adjust the output to accommodate for this excess if you choose to proceed with the current parameters,
      as downstream functions such as `compute_variant_sim`, `compute_input_score`, and `create_subsets` all rely on `n_variants`.")
    } else {
      keep.nonsig = sample(df_gwas.filtered.labeled %>% dplyr::filter(label == "nonsignificant") %>% dplyr::pull(SNP), size = n_to_sample)
      df_gwas.filtered.labeled = df_gwas.filtered.labeled %>%
        dplyr::mutate(included_in_study = dplyr::if_else(label == "nonsignificant" & !(SNP %in% keep.nonsig),
                                           FALSE,
                                           included_in_study))
    }
  }
  return(df_gwas.filtered.labeled)
}
