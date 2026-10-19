# Create subsets for DisCO-VG
# Adapted from Andre's code

# Elysia Chou, in collaboration with Andre Guerra

####### Generating Optimized Subset List #########
# Overview
# 1. Define sampling parameters
# 2. Generate optimized and balanced subsets
# 3. Save list

#' Define sampling parameters
#'
#' @param sample_pool_VG A list of data frames containing variant-gene pairs for each label (positive, prox_nonsig, nonsig)
#' @param MVG The full multi-variant-gene matrix (input scores)
#' @param n_repred How many times to sample each VG-pair
#' @param pos_vg_per_run Positive VG pairs per run
#' @param prox_ns_vg_per_run Proximal non-sig pairs per run
#' @param ns_vg_per_run Non-sig pairs per run
#' @param mat_size Size of mini-MVG matrix
#' @param sig_target_pct proportion of significant variants to include in mini matrix
#' @param prox_ns_target_pct proportion of proximal non-significant variants to include in mini matrix
#' @param seed Random seed for reproducibility
#' @return List of relevant sampling parameters to be used in subset generation
#' @export
list_mini_matrix <- function(sample_pool_VG, MVG,
                             n_repred = 10,
                             pos_vg_per_run = 40,
                             prox_ns_vg_per_run = 10,
                             ns_vg_per_run = 10,
                             mat_size = 1200,
                             sig_target_pct = 0.2,
                             prox_ns_target_pct = 0.1,
                             seed = 888) {
  # Calculate optimal number of runs based on significant pairs
  n_runs = ceiling(n_repred * nrow(sample_pool_VG[["positive"]]) / pos_vg_per_run)
  message("Will generate ", n_runs, " subsets")

  # Sample from these categories and keep track of original indices
  set.seed(seed)
  prox_ns_sample_idx = sample(1:nrow(sample_pool_VG[["prox_nonsig"]]), round(n_runs/n_repred*prox_ns_vg_per_run))
  ns_sample_idx = sample(1:nrow(sample_pool_VG[["nonsig"]]), round(n_runs/n_repred*ns_vg_per_run))

  # Number of VG pairs to sample per label
  prox_ns_vg = sample_pool_VG[["prox_nonsig"]][prox_ns_sample_idx, ]
  # prox_ns_vg$orig_idx = prox_ns_sample_idx
  ns_vg = sample_pool_VG[["nonsig"]][ns_sample_idx, ]
  # pos_vg = sample_pool_VG[["positive"]]

  # NOTE: Unnecessary.
  # # Initialize tracking arrays
  # pos_subsets = matrix(0, nrow = n_runs, ncol = pos_vg_per_run)
  # prox_ns_subsets = matrix(0, nrow = n_runs, ncol = prox_ns_vg_per_run)
  # ns_subsets = matrix(0, nrow = n_runs, ncol = ns_vg_per_run)

  # Sample for each class
  pos_result = sample_with_balance(seq_along(sample_pool_VG[["positive"]]$variant), pos_vg_per_run, n_runs, n_repred)
  prox_ns_result = sample_with_balance(prox_ns_sample_idx, prox_ns_vg_per_run, n_runs, n_repred)
  ns_result = sample_with_balance(ns_sample_idx, ns_vg_per_run, n_runs, n_repred)

  pos_subsets = pos_result$subset_matrix
  prox_ns_subsets = prox_ns_result$subset_matrix
  ns_subsets = ns_result$subset_matrix

  # Update counts in main index
  # sample_pool_VG[["positive"]]$count = 0
  sample_pool_VG[["prox_nonsig"]]$count = 0
  sample_pool_VG[["nonsig"]]$count = 0

  sample_pool_VG[["positive"]]$count = pos_result$actual_counts
  sample_pool_VG[["prox_nonsig"]]$count[prox_ns_sample_idx] = prox_ns_result$actual_counts
  sample_pool_VG[["nonsig"]]$count[ns_sample_idx] = ns_result$actual_counts

  # Verify balance
  message("Distribution of sampling counts:")
  table_counts = table(sample_pool_VG[["positive"]]$count)
  message("Significant pairs: ", paste(names(table_counts), table_counts, collapse=", "))

  message("Generating ", n_runs, " subsets...")
  subsets = vector(mode = "list", length = n_runs)

  for (i in 1:n_runs) {
    # Get indices for this run
    pos_idx = pos_subsets[i, pos_subsets[i,] > 0]
    prox_ns_idx = prox_ns_sample_idx[prox_ns_subsets[i, prox_ns_subsets[i,] > 0]]
    ns_idx = ns_sample_idx[ns_subsets[i, ns_subsets[i,] > 0]]

    # Get variants and genes for these VG pairs
    run_variants = c(
      sample_pool_VG[["positive"]]$variant[pos_idx],
      sample_pool_VG[["prox_nonsig"]]$variant[prox_ns_idx],
      sample_pool_VG[["nonsig"]]$variant[ns_idx]
    )
    run_genes = c(
      sample_pool_VG[["positive"]]$gene[pos_idx],
      sample_pool_VG[["prox_nonsig"]]$gene[prox_ns_idx],
      sample_pool_VG[["nonsig"]]$gene[ns_idx]
    )

    run_variants = unique(run_variants)
    run_genes = unique(run_genes)

    # Calculate how many additional variants we need of each class
    sig_vars = unique(sample_pool_VG[["positive"]]$variant)
    prox_ns_vars = unique(prox_ns_vg$variant)
    ns_vars = unique(ns_vg$variant)

    # Count how many of each we already have
    sig_in_run = sum(run_variants %in% sig_vars)
    prox_ns_in_run = sum(run_variants %in% prox_ns_vars)

    # Calculate how many more we need
    target_sig = round(mat_size * sig_target_pct)
    target_prox_ns = round(mat_size * prox_ns_target_pct)

    # Sample additional variants if needed
    set.seed(seed + i)  # Change seed for each run to ensure variability
    additional_sig = character(0)
    if (sig_in_run < target_sig) {
      additional_sig = sample(sig_vars[!sig_vars %in% run_variants],
                              size = min(target_sig - sig_in_run,
                                        sum(!sig_vars %in% run_variants)))
    }

    additional_prox_ns = character(0)
    if (prox_ns_in_run < target_prox_ns) {
      additional_prox_ns = sample(prox_ns_vars[!prox_ns_vars %in% run_variants],
                                  size = min(target_prox_ns - prox_ns_in_run,
                                            sum(!prox_ns_vars %in% run_variants)))
    }

    # Add non-sig variants to fill the matrix
    total_vars_so_far = length(run_variants) + length(additional_sig) + length(additional_prox_ns)
    n_remaining_vars = mat_size - total_vars_so_far

    additional_ns = character(0)
    if (n_remaining_vars > 0) {
      potential_ns = setdiff(rownames(MVG), c(run_variants, additional_sig, additional_prox_ns))
      additional_ns = sample(potential_ns, size = min(n_remaining_vars, length(potential_ns)))
    }

    # Combine all variants
    final_variants = c(run_variants, additional_sig, additional_prox_ns, additional_ns)

    # Add genes to fill the matrix if needed
    n_genes_needed = mat_size - length(run_genes)
    additional_genes = character(0)
    if (n_genes_needed > 0) {
      additional_genes = sample(setdiff(colnames(MVG), run_genes),
                                size = min(n_genes_needed, length(setdiff(colnames(MVG), run_genes))))
    }
    final_genes = c(run_genes, additional_genes)

    # Create set_to_zero data frame
    set_to_zero = data.frame(
      variant = c(
        sample_pool_VG[["positive"]]$variant[pos_idx],
        sample_pool_VG[["prox_nonsig"]]$variant[prox_ns_idx],
        sample_pool_VG[["nonsig"]]$variant[ns_idx]
      ),
      gene = c(
        sample_pool_VG[["positive"]]$gene[pos_idx],
        sample_pool_VG[["prox_nonsig"]]$gene[prox_ns_idx],
        sample_pool_VG[["nonsig"]]$gene[ns_idx]
      )
    )

    # Store in subset list
    subsets[[i]] = list(
      variants = final_variants[1:mat_size],  # Ensure consistent size
      genes = final_genes[1:mat_size],        # Ensure consistent size
      set_to_zero = set_to_zero,
      subset_ID = i
    )

    # Progress indicator
    if (i %% 100 == 0) message("Completed ", i, " of ", n_runs, " subsets")
  }
  return(subsets)
}

# Helper function to create balanced samples
sample_with_balance <- function(indices, pairs_per_run, n_runs, n_repred, seed = 999) {
  set.seed(seed)
  n_pairs = length(indices)
  count_vec = rep(0, n_pairs)

  # Pre-allocate result matrix
  result = matrix(0, nrow = n_runs, ncol = pairs_per_run)

  # First pass: distribute using weighted sampling
  for (run in 1:n_runs) {
    # Find pairs that haven't reached max count
    eligible = which(count_vec < n_repred)

    # If we've sampled everything n_repred times, break
    if (length(eligible) == 0) break

    # Calculate weights inversely proportional to current counts
    # Adding a small constant (0.1) to ensure non-zero probabilities
    p_weights = n_repred - count_vec[eligible] + 0.1

    # Normalize weights to probabilities
    probs = p_weights / sum(p_weights)

    # Sample without replacement using weighted probability
    n_to_sample = min(pairs_per_run, length(eligible))
    selected_indices = sample(eligible, size = n_to_sample,
                               prob = probs, replace = FALSE)

    # Pad with zeros if needed
    if (length(selected_indices) < pairs_per_run) {
      selected_indices = c(selected_indices, rep(0, pairs_per_run - length(selected_indices)))
    }

    result[run, ] = selected_indices
    count_vec[selected_indices[selected_indices > 0]] =
      count_vec[selected_indices[selected_indices > 0]] + 1
  }

  # Second pass: fill any remaining slots with weighted sampling
  remaining_slots = which(result == 0, arr.ind = TRUE)
  if (nrow(remaining_slots) > 0) {
    for (i in 1:nrow(remaining_slots)) {
      eligible = which(count_vec < n_repred)
      if (length(eligible) > 0) {
        # Use weighted sampling for remaining slots too
        p_weights = n_repred - count_vec[eligible] + 0.1
        probs = p_weights / sum(p_weights)
        selected_idx = sample(eligible, 1, prob = probs)

        result[remaining_slots[i,1], remaining_slots[i,2]] = selected_idx
        count_vec[selected_idx] = count_vec[selected_idx] + 1
      }
    }
  }

  # Final check - identify any pairs that haven't been sampled enough
  under_sampled = which(count_vec < n_repred)
  if (length(under_sampled) > 0) {
    message("Note: ", length(under_sampled), " pairs were sampled fewer than ",
            n_repred, " times")

    # Report distribution of sampling counts
    count_table = table(count_vec)
    message("Distribution of sampling counts: ",
            paste(names(count_table), count_table, collapse=", "))
  }

  return(list(subset_matrix = result, actual_counts = count_vec))
}
