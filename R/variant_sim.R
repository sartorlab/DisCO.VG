#' @title Variant-Variant Similarity Matrix (MVV) from ATAC-Seq Signal
#'
#' @description
#' Functions for constructing a variant-variant similarity matrix (MVV)
#' based on ATAC-Seq peak signal correlation across samples.
#'
#' @importFrom utils head
#'
#' @name variant_sim
NULL

# ===========================================================================
# create_mvv_from_peak_data
# ===========================================================================

#' Create an MVV directly from ATAC peak data
#'
#' Maps each variant to its overlapping ATAC peak, constructs a
#' variant x sample matrix, and computes the full MVV in a single pass.
#'
#' @param atac_matrix A numeric matrix of ATAC-Seq peak signal. Rows are
#'   peak IDs (e.g., \code{"chr1_100_600"}), columns are samples.
#' @param variants A data.frame with at least columns \code{SNP} and
#'   \code{peak_ID}. \code{peak_ID} must match row names of
#'   \code{atac_matrix}.
#' @param method Character. Correlation method passed to \code{\link[stats]{cor}}.
#'   Default \code{"pearson"}.
#' @param output_path Optional file path to save the MVV as an RDS file.
#'   Default \code{NULL} (no file saved).
#'
#' @return A symmetric numeric matrix (n_variants x n_variants) of pairwise
#'   correlations. Row and column names are variant SNP IDs.
#'
#' @details
#' Variants whose \code{peak_ID} is not found in \code{atac_matrix} are
#' filtered out with a warning. Variants mapping to the same peak will have
#' identical signal vectors and thus a correlation of 1.0.
#'
#' @examples
#' \dontrun{
#' atac <- as.matrix(read.csv("atac_peak_matrix.csv", row.names = 1))
#' variants <- read.csv("filtered_variants.csv")
#' mvv <- create_mvv_from_peak_data(atac, variants, output_path = "mvv.rds")
#' }
#'
#' @export
create_mvv_from_peak_data <- function(atac_matrix, variants, method = "pearson",
                                      output_path = NULL) {

  stopifnot(
    is.matrix(atac_matrix),
    is.numeric(atac_matrix),
    is.data.frame(variants),
    all(c("SNP", "peak_ID") %in% colnames(variants))
  )

  # Filter out variants with peak_IDs not in atac_matrix
  missing_peaks = setdiff(variants$peak_ID, rownames(atac_matrix))
  if (length(missing_peaks) > 0) {
    warning(
      "Filtering out ", sum(variants$peak_ID %in% missing_peaks),
      " variant(s) with peak_ID(s) not found in atac_matrix: ",
      paste(head(missing_peaks, 5), collapse = ", "),
      if (length(missing_peaks) > 5) ", ..." else ""
    )
    variants = variants[!variants$peak_ID %in% missing_peaks, , drop = FALSE]
  }

  stopifnot(nrow(variants) >= 2)
  variant_sample_matrix = atac_matrix[variants$peak_ID, , drop = FALSE]
  rownames(variant_sample_matrix) = variants$SNP

  mvv = compute_MVV(variant_sample_matrix, method = method, output_path = output_path)

  return(mvv)
}


# ===========================================================================
# create_mvv_from_MVV
# ===========================================================================

#' Create an MVV in parallel from chunked subsets
#'
#' Splits variants into chunks and computes within-chunk and between-chunk
#' correlation blocks in parallel using \code{\link[parallel]{mclapply}},
#' then assembles the full MVV.
#'
#' @param atac_matrix A numeric matrix of ATAC-Seq peak signal. Rows are
#'   peak IDs, columns are samples.
#' @param variants A data.frame with at least columns \code{SNP} and
#'   \code{peak_ID}.
#' @param n_chunks Integer. Number of chunks to split variants into.
#'   Default \code{4}.
#' @param n_cores Integer. Number of cores for parallel execution.
#'   Default \code{1} (sequential).
#' @param method Character. Correlation method. Default \code{"pearson"}.
#' @param output_path Optional file path to save the MVV as an RDS file.
#'   Default \code{NULL} (no file saved).
#'
#' @return A symmetric numeric matrix (n_variants x n_variants) of pairwise
#'   correlations.
#'
#' @details
#' Variants whose \code{peak_ID} is not found in \code{atac_matrix} are
#' filtered out with a warning.
#'
#' Both within-chunk and between-chunk correlations are computed, so the
#' result is identical to \code{create_mvv_from_peak_data}. The chunked
#' approach reduces per-worker memory and enables parallel speedup.
#'
#' Parallelization is OS-aware: uses \code{\link[parallel]{mclapply}}
#' (fork-based) on Unix/macOS and \code{\link[parallel]{parLapply}}
#' (PSOCK cluster) on Windows. Falls back to sequential \code{lapply}
#' when \code{n_cores = 1}.
#'
#' Falls back to \code{create_mvv_from_peak_data} when \code{n_chunks = 1}
#' or the variant set is small.
#'
#' @examples
#' \dontrun{
#' atac <- as.matrix(read.csv("atac_peak_matrix.csv", row.names = 1))
#' variants <- read.csv("filtered_variants.csv")
#' mvv <- create_mvv_from_MVV(atac, variants, n_chunks = 4, n_cores = 4)
#' }
#'
#' @export
create_mvv_from_MVV <- function(atac_matrix,
                                variants,
                                n_chunks = 4L,
                                n_cores = 1L,
                                method = "pearson",
                                output_path = NULL) {

  stopifnot(
    is.matrix(atac_matrix),
    is.numeric(atac_matrix),
    is.data.frame(variants),
    all(c("SNP", "peak_ID") %in% colnames(variants)),
    n_chunks >= 1L,
    n_cores >= 1L
  )

  # Filter out variants with peak_IDs not in atac_matrix
  missing_peaks = setdiff(variants$peak_ID, rownames(atac_matrix))
  if (length(missing_peaks) > 0) {
    warning(
      "Filtering out ", sum(variants$peak_ID %in% missing_peaks),
      " variant(s) with peak_ID(s) not found in atac_matrix: ",
      paste(head(missing_peaks, 5), collapse = ", "),
      if (length(missing_peaks) > 5) ", ..." else ""
    )
    variants = variants[!variants$peak_ID %in% missing_peaks, , drop = FALSE]
  }

  stopifnot(nrow(variants) >= 2)

  n_variants = nrow(variants)

  # Fall back to single-pass for trivial cases
  if (n_chunks <= 1L || n_variants < 10) {
    message("Using single-pass compute_MVV.")
    return(create_mvv_from_peak_data(atac_matrix, variants, method = method,
                                     output_path = output_path))
  }

  # Build full variant x sample matrix
  variant_sample_matrix = atac_matrix[variants$peak_ID, , drop = FALSE]
  rownames(variant_sample_matrix) = variants$SNP
  all_snps = variants$SNP

  # Split variant indices into chunks
  chunk_indices = split(
    seq_len(n_variants),
    cut(seq_len(n_variants), n_chunks, labels = FALSE)
  )

  message(
    "Computing MVV: ", n_variants, " variants, ",
    n_chunks, " chunks, ", n_cores, " core(s)"
  )

  # Enumerate all chunk pairs (upper triangle + diagonal)
  chunk_pairs = list()
  for (i in seq_along(chunk_indices)) {
    for (j in i:length(chunk_indices)) {
      chunk_pairs = c(chunk_pairs, list(c(i, j)))
    }
  }

  # Function to compute one correlation block
  compute_block = function(pair) {
    ci = pair[1]
    cj = pair[2]

    if (ci == cj) {
      block = stats::cor(
        t(variant_sample_matrix[chunk_indices[[ci]], , drop = FALSE]),
        method = method
      )^2
    } else {
      block = stats::cor(
        t(variant_sample_matrix[chunk_indices[[ci]], , drop = FALSE]),
        t(variant_sample_matrix[chunk_indices[[cj]], , drop = FALSE]),
        method = method
      )^2
    }

    list(i = ci, j = cj, block = block)
  }

  # Compute correlation blocks in parallel (auto-detect OS)
  is_windows = .Platform$OS.type == "windows"

  if (n_cores > 1 && is_windows) {
    # Windows: use PSOCK cluster with parLapply
    cl = parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Export required objects to workers
    parallel::clusterExport(cl, c("variant_sample_matrix", "chunk_indices", "method"),
                            envir = environment())

    block_results = parallel::parLapply(cl, chunk_pairs, compute_block)

  } else if (n_cores > 1) {
    # Unix/macOS: use mclapply (fork-based)
    block_results = parallel::mclapply(chunk_pairs, compute_block, mc.cores = n_cores)

  } else {
    # Sequential fallback
    block_results = lapply(chunk_pairs, compute_block)
  }

  # Assemble full MVV from blocks
  mvv = matrix(NA_real_, nrow = n_variants, ncol = n_variants,
               dimnames = list(all_snps, all_snps))

  for (res in block_results) {
    idx_i = chunk_indices[[res$i]]
    idx_j = chunk_indices[[res$j]]
    mvv[idx_i, idx_j] = res$block
    if (res$i != res$j) {
      mvv[idx_j, idx_i] = t(res$block)
    }
  }

  if (!is.null(output_path)) {
    saveRDS(mvv, output_path)
    message("MVV saved to: ", output_path)
  }

  return(mvv)
}


# ===========================================================================
# compute_MVV
# ===========================================================================

#' Compute a variant-variant similarity matrix
#'
#' Core function. Takes a variant x sample matrix and computes pairwise
#' Pearson correlation across all variants.
#'
#' @param variant_sample_matrix A numeric matrix with variants as rows and
#'   samples as columns. Row names should be variant IDs.
#' @param method Character. Correlation method passed to \code{\link[stats]{cor}}.
#'   Default \code{"pearson"}.
#'
#' @return A symmetric numeric matrix (n_variants x n_variants) with
#'   correlation values. Row and column names are variant IDs.
#'
#' @param output_path Optional file path to save the MVV as an RDS file.
#'   Default \code{NULL} (no file saved).
#'
#' @examples
#' mat <- matrix(rnorm(20), nrow = 4, dimnames = list(
#'   paste0("rs", 1:4), paste0("sample_", 1:5)
#' ))
#' mvv <- compute_MVV(mat)
#'
#' @export
compute_MVV <- function(variant_sample_matrix, method = "pearson", output_path = NULL) {

  stopifnot(
    is.matrix(variant_sample_matrix),
    is.numeric(variant_sample_matrix),
    nrow(variant_sample_matrix) >= 2,
    ncol(variant_sample_matrix) >= 2
  )

  # cor() operates on columns; transpose so variants are correlated
  mvv = stats::cor(t(variant_sample_matrix), method = method)^2

  if (!is.null(output_path)) {
    saveRDS(mvv, output_path)
    message("MVV saved to: ", output_path)
  }

  return(mvv)
}
