#' @title Run CMMC algorithm to generate raw DisCO-VG scores
#'
#' @description
#' Given the main input score matrix and the side matrices providing supporting
#' info, compute raw DisCO-VG scores using the CMMC algorithm on multiple bagged
#' subsets of the main data. Merges the results of these many runs at the end.
#'
#' @importFrom dplyr filter
#' @importFrom data.table rbindlist
#'
#' @name run_CMMC
NULL

# Create mini-matrices and run CMMC algorithm for each mini-matrix.
# Merge all raw results at the end.

# Elysia Chou, in collaboration with Andre Guerra

#' Create mini-matrices for a given variant-gene subset
#'
#' Create mini-matrices for a given subset of variants and genes to pass on to CMMC algorithm.
#'
#' @param mini_matrix_info A list with the following keys:
#'  - subset_ID
#'  - variants (which variants are included in the subset),
#'  - genes (which genes are included in the subset),
#'  - set_to_zero (dataframe with columns variant and gene that specifies which pairs to impute in CMMC algorithm).
#' @param Matrix_info A list with the following keys:
#'  - MVG (the full variant-gene matrix),
#'  - MVV (the full variant-variant matrix),
#'  - MGG (the full gene-gene matrix),
#'  - peak_data (boolean specifying whether the MVV is peak x sample data or a precomputed similarity matrix).
#' @param full_variants A list of the GWAS variants being analyzed,
#'  containing the information important for constructing the variant-variant
#'  similarity matrix (mvv) from peak data.
#' @param gene_list List of genes being considered in the DisCO-VG study.
#'  Defaults to built-in gene list that is filtered for genes that have at least one GO
#'  BP (Biological Process) term.
#' @return A list with the following keys:
#'  - mvg (the mini variant-gene matrix),
#'  - mvv (the mini variant-variant matrix),
#'  - mgg (the mini gene-gene matrix).
#'
#' @export
create_mini_matrix <- function(mini_matrix_info, Matrix_info,
                               full_variants, gene_list = NULL) {

  # Load default gene_list if necessary
  if (is.null(gene_list)) {
    gene_list = load_gene_list()
  }
  # TODO: catch errors when not all variants or genes in subset are in the full matrix.
  # TODO: see if you still need full_variants if peak_data is FALSE.

  # --- Create subsets (lowercase = subset) ---
  # -- Main matrix --
  # Proposed code assuming subsets is a list of lists with variants and genes as keys
  mvg = matrix(0,
                nrow = length(mini_matrix_info$variants),
                ncol = length(mini_matrix_info$genes),
                dimnames = list(mini_matrix_info$variants, mini_matrix_info$genes))

  for (variant in mini_matrix_info$variants) {
      mvg[variant,mini_matrix_info$genes] = Matrix_info$MVG[variant,mini_matrix_info$genes]
  }

  # -- variant-variant matrix --
  # Specify that default is TRUE
  # Order here is very important
  df_vars = full_variants[match(mini_matrix_info$variants, full_variants$SNP),]

  if (Matrix_info$peak_data) {
      # TODO: add option for similarity method; using default for now.
      mvv = create_mvv_from_peak_data(atac_matrix = Matrix_info$MVV,
                                       variants = df_vars)
  } else {
      mvv = create_mvv_from_MVV(atac_matrix = Matrix_info$MVV,
                                 variants = df_vars)
  }

  # -- gene-gene matrix --
  # GOBP
  mgg = create_mgg(mini_matrix_info, Matrix_info$MGG, gene_list)

  return(list(mvg = mvg, mvv = mvv, mgg = mgg))
}

#' Compute raw DisCO-VG scores for a given mini-matrix
#'
#' Wrapper function to obtain a dataframe of DisCO-VG input and raw ouptut scores
#' for a given mini-matrix. This function will first run the CMMC algorithm and
#' compute raw DisCO-VG scores for the specified variant-gene pairs in
#' `mini_matrix_info[["set_to_zero"]]`. (TODO: PROVIDE MORE DETAILS HERE ABOUT WHY WE DO THIS)
#'
#' @param mini_matrix_info A list with the following keys:
#'  - subset_ID
#'  - variants (which variants are included in the subset),
#'  - genes (which genes are included in the subset),
#'  - set_to_zero (dataframe with columns variant and gene that specifies which pairs to impute in CMMC algorithm).
#' @param mini_matrices A list with the following keys:
#'  - mvg (the mini variant-gene matrix),
#'  - mvv (the mini variant-variant matrix),
#'  - mgg (the mini gene-gene matrix).
#' @param CMMC_config A list with the following keys:
#'  - x_epsilon
#'  - x_scale
#'  - y_epsilon
#'  - y_scale
#'  - iter_num
#' @param include_outputs List of Booleans about what CMMC output parameters to
#'  include for the helper function `list_outputs()`:
#'  - `xy_hat`: Boolean to indicate whether to include the CMMC-completed matrix
#'    in the list of outputs.
#'  - `learned_weights`: Boolean to indicate whether to include the per-slice weights
#'    learned during optimization (`w_x`, `w_y`) in the list of outputs. Only applicable
#'    to CTMC applications.
#'  - `num_slices`: Boolean to indicate whether to include the number of variant-variant
#'    and gene-gene similarity slices (`nU`, `nZ`) in the list of outputs. Only applicable
#'    to CTMC applications, as CTMC uses similarity tensors (i.e., multiple slices
#'    of similarity matrices).
#'  - `iters_used`: Boolean to indicate whether to include the number of iterations
#'    used during the CMMC optimization in the list of outputs.
#'  - `F_final`: Boolean to indicate whether to include the final objective function
#'    value in the list of outputs.
#'  - `F_trace`: Boolean to indicate whether to include the objective function
#'    trace in the list of outputs.
#'  - `converged`: Boolean indicating whether the CMMC algorithm converged or not. If FALSE, more diagnostic
#'    outputs are included in the final list.
#' @return A list of outputs from the CMMC algorithm; see the `list_outputs()` function.
#' This function will by default return:
#'  - `df_scores`: A dataframe with raw DisCO-VG scores for the specified
#'    variant-gene pairs in `mini_matrix_info[["set_to_zero"]]`.
#'    The dataframe has the following columns:
#'    - index (the subset ID)
#'    - variant
#'    - gene
#'    - input_score (the original score in the mini variant-gene matrix for the specified variant-gene pair)
#'    - output_score (the raw, imputed score from the CMMC output matrix for the specified variant-gene pair)
#'  - `F_trace`: A vector tracing the objective function for running CMMC on this mini-matrix.
#'  - `converged`: whether the CMMC algorithm converged on this mini-matrix or not.
#'    If FALSE, more diagnostic outputs are included in the final list.
#'
#' @export
compute_raw_score <- function(mini_matrix_info, mini_matrices,
                              CMMC_config = list(
                                x_epsilon = 1.0,
                                x_scale = 1.0,
                                y_epsilon = 1.0,
                                y_scale = 1.0,
                                iter_num = 2
                              ),
                              include_outputs = list(
                                xy_hat = FALSE,
                                learned_weights = FALSE,
                                num_slices = FALSE,
                                iters_used = FALSE,
                                F_final = FALSE,
                                F_trace = TRUE,
                                converged = TRUE
                              )) {
  # Run CMMC algorithm
  CMMC_output = run_CMMC(mini_matrix_info, mini_matrices, CMMC_config)

  # Organize scores into a list
  outputs = list_outputs(CMMC_output, mini_matrix_info, mini_matrices,
                         xy_hat = include_outputs$xy_hat,
                         learned_weights = include_outputs$learned_weights,
                         num_slices = include_outputs$num_slices,
                         iters_used = include_outputs$iters_used,
                         F_final = include_outputs$F_final,
                         F_trace = include_outputs$F_trace,
                         converged = include_outputs$converged)
  return(outputs)
}

#' Merge raw DisCO-VG scores from multiple subsets into a single dataframe
#'
#' If you compute multiple raw DisCO-VG output scores for a given variant-gene pair
#' by running the CMMC algorithm on numerous subsets (mini-matrices),
#' this function merges a list of dataframes returned by `compute_raw_score` into a single dataframe.
#' If any run on a given subset returned an error instead of a dataframe, the function will flag these runs for further review.
#' @param output_list A list of lists of dataframes returned by `compute_raw_score` for multiple runs of the CMMC algorithm
#'   on different subsets (mini-matrices).
#' @return A list with the following keys:
#'  - final_df: a single dataframe with raw DisCO-VG scores for all subsets.
#'  - errors: a list of subsets that returned errors instead of dataframes, if any; otherwise NULL.
#'
#' @export
merge_raw_scores <- function(output_list) {
  indices_failed = which(sapply(seq_along(output_list), function(x) {return(class(output_list[[x]]$df_scores) != "data.frame")}))
  if (length(indices_failed) > 0) {
    warning(paste("Errors found in the following indices; see 'errors' key for more details:\n",
            paste(indices_failed, collapse = ", ")))
    final_df = data.table::rbindlist(lapply(output_list[-indices_failed], function(x) {return(x$df_scores)}))
    # Save errors for debugging
    errors = lapply(output_list[indices_failed], function(x) {return(x$df_scores)})
    return(list(final_df = final_df, errors = errors))
  } else {
    final_df = data.table::rbindlist(lapply(output_list, function(x) {return(x$df_scores)}))
    return(list(final_df = final_df, errors = NULL))
  }
}

#' Merge objective function traces from multiple subsets into a single dataframe
#'
#' @param output_list A list of lists of outputs returned by `compute_raw_score`
#'  for multiple runs of the CMMC algorithm on different subsets (mini-matrices).
#'  In particular, the sublists must have `F_trace` as a key.
#' @return A list with the following keys:
#'  - `df.obj_fxn`: dataframe with all the objective function values
#'  -  errors: a list of subsets that returned errors or do not have `F_trace` as a key.
#' @export
merge_F_trace <- function(output_list) {
  indices_failed = which(sapply(seq_along(output_list), function(x) {return(!("F_trace" %in% names(output_list[[x]])))}))
  if (length(indices_failed) > 0) {
    warning(paste("No objective function trace recorded in the following subsets; see 'errors' key for the subset IDs:\n",
                  paste(indices_failed, collapse = ", ")))
    vals = unlist(lapply(output_list[-indices_failed], function(x) {x$F_trace}))
    ids = rep(x = unlist(lapply(output_list[-indices_failed], function(x) {x$df_scores$index[1]})),
              times = unlist(lapply(output_list[-indices_failed], function(x) {length(x$F_trace)})))
    pos = unlist(lapply(output_list[-indices_failed], function(x) {seq_along(x$F_trace)}))

    df.obj_fxn = data.frame(
      F_trace = vals,
      subset_ID = ids,
      iteration = pos
    )

    errors = names(output_list)[indices_failed]
    return(list(df.obj_fxn = df.obj_fxn, errors = errors))
  } else {
    vals = unlist(lapply(output_list, function(x) {x$F_trace}))
    ids = rep(x = unlist(lapply(output_list, function(x) {x$df_scores$index[1]})),
              times = unlist(lapply(output_list, function(x) {length(x$F_trace)})))
    pos = unlist(lapply(output_list, function(x) {seq_along(x$F_trace)}))

    df.obj_fxn = data.frame(
      F_trace = vals,
      subset_ID = ids,
      iteration = pos
    )
    return(list(df.obj_fxn = df.obj_fxn, errors = NULL))
  }
}

#' Run the CMMC algorithm to compute raw scores for a given mini-matrix
#'
#' @param mini_matrix_info A list with the following keys:
#'  - subset_ID
#'  - variants (which variants are included in the subset),
#'  - genes (which genes are included in the subset),
#'  - set_to_zero (dataframe with columns variant and gene that specifies which pairs to impute in CMMC algorithm).
#' @param mini_matrices A list with the following keys:
#'  - mvg (the mini variant-gene matrix),
#'  - mvv (the mini variant-variant matrix),
#'  - mgg (the mini gene-gene matrix).
#' @param CMMC_config A list with the following keys:
#'  - x_epsilon
#'  - x_scale
#'  - y_epsilon
#'  - y_scale
#'  - iter_num
#' @return A list of outputs from the `CMMC()` function:
#'  - `xy_hat`: Completed mini-matrix as a result of the CMMC algorithm.
#'  - `w_x`, `w_y`: In the CTMC extension, which uses similarity tensors representing
#'    slices of similarity matrices, `w_x` and `w_y` are the per-slice weights learned
#'    during optimization.
#'  - `nU`, `nZ`: in the CTMC extension, which uses similarity tensors, `nU` and `nZ`
#'    represent the number of variant-variant similarity slices and gene-gene
#'    similarity slices, respectively.
#'  - `iters_used`: number of iterations until convergence, or max number of iterations.
#'  - `F_final`: final value of the objective function
#'  - `F_trace`: vector that tracks the objective function over `iters_used` iterations.
#'
#' @export
run_CMMC <- function(mini_matrix_info, mini_matrices,
                     CMMC_config = list(
                       x_epsilon = 1.0,
                       x_scale = 1.0,
                       y_epsilon = 1.0,
                       y_scale = 1.0,
                       iter_num = 2
                     )) {
  # -- Record the score and set specific pairs in main matrix to zero --
  for (i in 1:nrow(mini_matrix_info$set_to_zero)) {
    variant = mini_matrix_info$set_to_zero[i, "variant"]
    gene = mini_matrix_info$set_to_zero[i, "gene"]
    mini_matrices$mvg[variant, gene] = 0
  }

  # -- Run CMMC --
  cmmc_result = CMMC(mini_matrices$mvg,
                     mini_matrices$mvv,
                     mini_matrices$mgg,
                     CMMC_config$x_scale,
                     CMMC_config$y_scale,
                     CMMC_config$x_epsilon,
                     CMMC_config$y_epsilon,
                     CMMC_config$iter_num)

  # The output is a list with a matrix as one of its items;
  # the output matrix with no row or column names
  rownames(cmmc_result$xy_hat) = rownames(mini_matrices$mvg)
  colnames(cmmc_result$xy_hat) = colnames(mini_matrices$mvg)

  return(cmmc_result)
}

#' Extract raw DisCO-VG scores after running the CMMC algorithm for a given mini-matrix
#'
#' @param cmmc_result A list of outputs given by the CMMC algorithm, usually obtained
#'    by running `CMMC()` or the DisCO-VG-specific wrapper `run_CMMC()` on a given mini-matrix.
#' @param mini_matrix_info A list with the following keys:
#'  - subset_ID
#'  - variants (which variants are included in the subset),
#'  - genes (which genes are included in the subset),
#'  - set_to_zero (dataframe with columns variant and gene that specifies which pairs to impute in CMMC algorithm).
#' @param mini_matrices A list with the following keys:
#'  - mvg (the mini variant-gene matrix),
#'  - mvv (the mini variant-variant matrix),
#'  - mgg (the mini gene-gene matrix).
#' @return A dataframe with raw DisCO-VG scores for the specified variant-gene pairs in `mini_matrix_info[["set_to_zero"]]`. The dataframe has the following columns:
#'  - index (the subset ID)
#'  - variant
#'  - gene
#'  - input_score (the original score in the mini variant-gene matrix for the specified variant-gene pair)
#'  - output_score (the raw, imputed score from the CMMC output matrix for the specified variant-gene pair)
#'
#' @export
extract_raw_scores_from_CMMC <- function(cmmc_result,
                                         mini_matrix_info,
                                         mini_matrices) {
  # -- Create dataframe of original and computed DisCO-VG scores --
  df_scores = data.frame(index = integer(),
                         variant = character(),
                         gene = character(),
                         input_score = numeric(),
                         output_score = numeric(),
                         stringsAsFactors = FALSE)

  # -- Record the input and output scores --
  # for the variant-gene pairs that are re-imputed with CMMC
  for (i in 1:nrow(mini_matrix_info$set_to_zero)) {
    variant = mini_matrix_info$set_to_zero[i, "variant"]
    gene = mini_matrix_info$set_to_zero[i, "gene"]
    temp_df = data.frame(index = mini_matrix_info$subset_ID,
                         variant = variant,
                         gene = gene,
                         input_score = mini_matrices$mvg[variant, gene],
                         output_score = cmmc_result$xy_hat[variant, gene],
                         stringsAsFactors = FALSE)
    df_scores = rbind(df_scores, temp_df)
  }
  # return dataframe for future merging
  return(df_scores)
}

#' Create list of outputs after running CMMC on a given mini-matrix
#'
#' For each mini-matrix, decide what output parameters and scores you would like to save.
#' The most common example is saving the dataframe of raw DisCO-VG output scores
#' with their corresponding input as well as saving the objective function trace
#' `F_trace`.
#'
#' @param cmmc_result A list of outputs given by the CMMC algorithm, usually obtained
#'    by running `CMMC()` or the DisCO-VG-specific wrapper `run_CMMC()` on a given mini-matrix.
#' @param mini_matrix_info A list with the following keys:
#'  - subset_ID
#'  - variants (which variants are included in the subset),
#'  - genes (which genes are included in the subset),
#'  - set_to_zero (dataframe with columns variant and gene that specifies which pairs to impute in CMMC algorithm).
#' @param mini_matrices A list with the following keys:
#'  - mvg (the mini variant-gene matrix),
#'  - mvv (the mini variant-variant matrix),
#'  - mgg (the mini gene-gene matrix).
#' @param xy_hat Boolean to indicate whether to include the CMMC-completed matrix
#'  in the list of outputs.
#' @param learned_weights Boolean to indicate whether to include the per-slice weights
#'  learned during optimization (`w_x`, `w_y`) in the list of outputs. Only applicable
#'  to CTMC applications.
#' @param num_slices Boolean to indicate whether to include the number of variant-variant
#'  and gene-gene similarity slices (`nU`, `nZ`) in the list of outputs. Only applicable
#'  to CTMC applications, as CTMC uses similarity tensors (i.e., multiple slices
#'  of similarity matrices).
#' @param iters_used Boolean to indicate whether to include the number of iterations
#'  used during the CMMC optimization in the list of outputs.
#' @param F_final Boolean to indicate whether to include the final objective function
#'  value in the list of outputs.
#' @param F_trace Boolean to indicate whether to include the objective function
#'  trace in the list of outputs.
#' @return A list of selected outputs:
#'  - `df_scores`: A dataframe of DisCO-VG input and raw output scores computed
#'    for a given mini-matrix. An output of the funciton `extract_raw_scores_from_CMMC()`.
#'    NOTE: If the algorithm didn't converge, this dataframe will be empty.
#'  - Optional:
#'    - `xy_hat`: Completed mini-matrix as a result of the CMMC algorithm.
#'    - `w_x`, `w_y`: In the CTMC extension, which uses similarity tensors representing
#'      slices of similarity matrices, `w_x` and `w_y` are the per-slice weights learned
#'      during optimization.
#'    - `nU`, `nZ`: in the CTMC extension, which uses similarity tensors, `nU` and `nZ`
#'      represent the number of variant-variant similarity slices and gene-gene
#'      similarity slices, respectively.
#'    - `iters_used`: number of iterations until convergence, or max number of iterations.
#'    - `F_final`: final value of the objective function
#'    - `F_trace`: vector that tracks the objective function over `iters_used` iterations.
#'    - other options that will receive an explanation at a later time (TODO)
#'
#' @export
list_outputs <- function(cmmc_result,
                         mini_matrix_info,
                         mini_matrices,
                         xy_hat = FALSE,
                         learned_weights = FALSE,
                         num_slices = FALSE,
                         iters_used = FALSE,
                         F_final = FALSE,
                         F_trace = TRUE,
                         fit_status = FALSE,
                         flags = FALSE,
                         failure = FALSE,
                         converged = TRUE,
                         stop_reason = FALSE,
                         solver = FALSE,
                         solve_status = FALSE,
                         relres = FALSE,
                         eig_floored = FALSE,
                         cg_iters = FALSE) {
  outputs = list()

  # If CMMC/CTMC didn't converge on a particular matrix setup, return an empty dataframe
  if (cmmc_result[["converged"]]) {
    outputs[["df_scores"]] = extract_raw_scores_from_CMMC(cmmc_result, mini_matrix_info, mini_matrices)
  } else {
    outputs[["df_scores"]] = df_scores = data.frame(index = integer(),
                                                    variant = character(),
                                                    gene = character(),
                                                    input_score = numeric(),
                                                    output_score = numeric(),
                                                    stringsAsFactors = FALSE)
    # Also include diagnostic information
    fit_status = TRUE
    flags = TRUE
    failure = TRUE
    converged = TRUE
    stop_reason = TRUE
    solver = TRUE
    solve_status = TRUE
    relres = TRUE
    eig_floored = TRUE
  }


  if (xy_hat) {
    outputs[["xy_hat"]] = cmmc_result[["xy_hat"]]
  }

  if (learned_weights) {
    outputs[["w_x"]] = cmmc_result[["w_x"]]
    outputs[["w_y"]] = cmmc_result[["w_y"]]
  }

  if (num_slices) {
    outputs[["nU"]] = cmmc_result[["nU"]]
    outputs[["nZ"]] = cmmc_result[["nZ"]]
  }

  if (iters_used) {
    outputs[["iters_used"]] = cmmc_result[["iters_used"]]
  }

  if (F_final) {
    outputs[["F_final"]] = cmmc_result[["F_final"]]
  }

  if (F_trace) {
    outputs[["F_trace"]] = cmmc_result[["F_trace"]]
  }

  if (fit_status) {
    outputs[["fit_status"]] = cmmc_result[["fit_status"]]
  }

  if (flags) {
    outputs[["flags"]] = cmmc_result[["flags"]]
  }

  if (failure) {
    outputs[["failure"]] = cmmc_result[["failure"]]
  }

  if (converged) {
    outputs[["converged"]] = cmmc_result[["converged"]]
  }

  if (stop_reason) {
    outputs[["stop_reason"]] = cmmc_result[["stop_reason"]]
  }

  if (solver) {
    outputs[["solver"]] = cmmc_result[["solver"]]
  }

  if (solve_status) {
    outputs[["solve_status"]] = cmmc_result[["solve_status"]]
  }

  if (relres) {
    outputs[["relres"]] = cmmc_result[["relres"]]
  }

  if (eig_floored) {
    outputs[["eig_floored"]] = cmmc_result[["eig_floored"]]
  }

  if (cg_iters) {
    outputs[["cg_iters"]] = cmmc_result[["cg_iters"]]
  }

  return(outputs)
}

#' Helper function to make mgg mini matrix
#'
#' Create mini matrix for gene-gene similarity, termed mgg, from the full
#' gene-gene similarity matrix, MGG. Note that MGG stores Entrez ID, while MVG
#' uses Ensembl ID.
#'
#' @param mini_matrix_info A list with the following keys:
#'  - subset_ID
#'  - variants (which variants are included in the subset),
#'  - genes (which genes are included in the subset),
#'  - set_to_zero (dataframe with columns variant and gene that specifies which pairs to impute in CMMC algorithm).
#' @param MGG Full gene-gene similarity matrix. The `DisCO.VG` package has a default
#'  MGG you can load for this purpose.
#' @param gene_list List of genes being considered in the DisCO-VG study.
#'  Defaults to built-in gene list that is filtered for genes that have at least one GO
#'  BP (Biological Process) term.
create_mgg <- function(mini_matrix_info, MGG, gene_list) {
  gene_idx = match(mini_matrix_info$genes, gene_list$ensembl_id)

  if (sum(is.na(gene_idx)) > 0) {
    warning("There seems to be a mismatch between the subset's", mini_matrix_info$subset_ID, "gene list and the provided gene list.\nRemoving NA genes...")
    gene_idx = gene_idx[!is.na(gene_idx)]
  }

  genes = as.character(gene_list$ensembl_id[gene_idx])

  mgg = as.matrix(MGG[genes, genes])

  return(mgg)
}
