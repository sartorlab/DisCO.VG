#' @title Process raw DisCO-VG output scores
#'
#' @description
#' Aggregate DisCO-VG output scores based on the various times a variant-gene pair was imputed.
#' For each variant-gene pair, the DisCO-VG output score is a summary statistic
#' computed from its raw scores from multiple, bagged imputations. This module
#' also then computes all relevant DisCO-VG scores -- input, output, and average
#' between input and output.
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
#' @return dataframe of DisCO-VG input and output scores where each row is a unique variant-gene pair
#' @export
process_output_score <- function(df_raw, variants, output_type = c("max", "overallTop10perc")) {
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

  if (output_type == "max") {
    df_output = compute_max_score(df_raw)
  } else if (output_type == "overallTop10perc") {
    df_output = compute_Top10perc(df_raw)
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
            any(startsWith(colnames(df_output), "DisCOVGscore")))

  df_final = tidyr::pivot_wider(df_output,
                                names_from = scoreType,
                                names_prefix = "DisCO.VG.",
                                values_from = dplyr::starts_with("DisCOVGscore"))

  df_final$DisCO.VG.avg_in_out = (df_final$DisCO.VG.input + df_final$DisCO.VG.output) / 2

  return(dplyr::rename(df_final, GWAS_label = label))
}

#' Compute max score from raw score
#'
#' @param df_raw.labeled merged dataframe of repredicted scores with label for each variant
#' @return dataframe where each row is a unique variant-gene pair with 90th percentile output score
compute_max_score <- function(df_raw.labeled) {
  df = tidyr::pivot_longer(df_raw.labeled,
                           cols = c("input_score","output_score"),
                           names_to = "scoreType",
                           names_pattern = "^(.*)_score$",
                           values_to = "DisCOVGscore")

  df = dplyr::summarize(df,
                        DisCOVGscore.max = max(DisCOVGscore),
                        .by = c(variant, label, gene, scoreType))

  return(df)
}

#' Compute 90th percentile score from raw score
#'
#' @param df_raw.labeled merged dataframe of repredicted scores with label for each variant
#' @return dataframe where each row is a unique variant-gene pair with 90th percentile output score
compute_Top10perc <- function(df_raw.labeled) {
  df = tidyr::pivot_longer(df_raw.labeled,
                           cols = c("input_score","output_score"),
                           names_to = "scoreType",
                           names_pattern = "^(.*)_score$",
                           values_to = "DisCOVGscore")

  df = dplyr::summarize(df,
                        DisCOVGscore.overallTop10perc = stats::quantile(DisCOVGscore, 0.9),
                        .by = c(variant, label, gene, scoreType))

  return(df)
}
