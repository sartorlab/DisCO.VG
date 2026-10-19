#' @title Process raw DisCO-VG output scores
#'
#' @description
#' Aggregate DisCO-VG output scores based on the various times a variant-gene pair was imputed.
#' For each variant-gene pair, the DisCO-VG score is a summary statistic
#' computed from either its raw input, output, and average between input and output scores
#' from multiple, bagged imputations. In the original study, avg(in, out) scores
#' are used as final DisCO-VG scores, reflecting data fusion from both input and output.
#'
#' @importFrom dplyr left_join join_by rename starts_with summarize
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom stats quantile
#'
#' @name compute_output_score
NULL

#' Process DisCO-VG output score from raw scores
#'
#' @param df_raw merged dataframe of repredicted scores (multiple scores per variant-gene pair
#' @param variants dataframe of variants used in model construction, with at minimum columns `SNP` and `label`
#' @param output_type specifies the summary statistic used as output score; currently only accepts max score or 90th percentile of scores for a given variant-gene pair.
#' @return dataframe of DisCO-VG input, output, and avg(in, out) scores where each group of 3 rows is a unique variant-gene pair
#' @export
process_output_score <- function(df_raw, variants, output_type = c("max", "90thPerc")) {
  # Validate output_type
  output_type = match.arg(output_type)

  # Validate necessary columns
  stopifnot(is.data.frame(variants),
            all(c("SNP", "label") %in% colnames(variants)),
            is.data.frame(df_raw),
            all(c("variant", "gene", "input_score", "output_score", "index") %in% colnames(df_raw)))

  df_raw = dplyr::left_join(df_raw,
                             variants[, c("SNP", "label")],
                             by = dplyr::join_by(variant == SNP))

  df_raw$avg_in_out_score = (df_raw$input_score + df_raw$output_score) / 2

  if (output_type == "max") {
    df_output = compute_max_score(df_raw)
  } else if (output_type == "90thPerc") {
    df_output = compute_90thPerc(df_raw)
  }

  return(df_output)
}

#' Provide finalized DisCO-VG scores
#'
#' Computes the recommended DisCO-VG avg(input, output) score after processing
#' output scores for each variant-gene pair.
#'
#' @param df_output dataframe of DisCO-VG input and output scores where each row is a unique variant-gene pair
#' @return dataframe with columns:
#'  - variant
#'  - gene
#'  - GWAS variant label
#'  - DisCO-VG input score
#'  - DisCO-VG output score
#'  - DisCO-VG average between input and output score
#' @export
DisCO.VG_scores <- function(df_output) {
  stopifnot(is.data.frame(df_output),
            all(c("variant", "gene", "label", "scoreType") %in% colnames(df_output)),
            any(startsWith(colnames(df_output), "DisCO.VG_score")))

  df_final = tidyr::pivot_wider(df_output,
                                names_from = scoreType,
                                names_prefix = "DisCO.VG.",
                                values_from = dplyr::starts_with("DisCO.VG_score"))

  # Order by DisCO.VG.avg_in_out
  df_final = df_final[order(df_final$DisCO.VG.avg_in_out, decreasing = TRUE),]
  return(dplyr::rename(df_final, GWAS_label = label))
}

#' Compute max score from raw score
#'
#' @param df_raw.labeled merged dataframe of repredicted scores with label for each variant
#' @return dataframe where each row is a unique variant-gene pair with 90th percentile output score
compute_max_score <- function(df_raw.labeled) {
  df = tidyr::pivot_longer(df_raw.labeled,
                           cols = c("input_score", "output_score", "avg_in_out_score"),
                           names_to = "scoreType",
                           names_pattern = "^(.*)_score$",
                           values_to = "DisCO.VG_score")

  df = dplyr::summarize(df,
                        DisCO.VG_score_max = max(DisCO.VG_score),
                        .by = c(variant, label, gene, scoreType))

  return(df)
}

#' Compute 90th percentile score from raw score
#'
#' @param df_raw.labeled merged dataframe of repredicted scores with label for each variant
#' @return dataframe where each row is a unique variant-gene pair with 90th percentile output score
compute_90thPerc <- function(df_raw.labeled) {
  df = tidyr::pivot_longer(df_raw.labeled,
                           cols = c("input_score", "output_score", "avg_in_out_score"),
                           names_to = "scoreType",
                           names_pattern = "^(.*)_score$",
                           values_to = "DisCO.VG_score")

  df = dplyr::summarize(df,
                        DisCO.VG_score_90thPerc = stats::quantile(DisCO.VG_score, 0.9),
                        .by = c(variant, label, gene, scoreType))

  return(df)
}
