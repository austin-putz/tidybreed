#' Select loci by count (n)
#'
#' @param genome Data frame of genome_meta table
#' @param n Number of SNPs to select
#' @param method Selection method: "random", "even", or "chromosome_even"
#' @return Logical vector indicating selected loci
#' @keywords internal
select_by_n <- function(genome, n, method) {

  n_loci <- nrow(genome)

  # Validate n
  if (n > n_loci) {
    stop(
      "n (", n, ") exceeds total loci (", n_loci, ")",
      call. = FALSE
    )
  }

  chip_vector <- rep(FALSE, n_loci)

  if (method == "random") {
    # Random selection
    selected_indices <- sample(1:n_loci, size = n, replace = FALSE)
    chip_vector[selected_indices] <- TRUE

  } else if (method == "even") {
    # Even spacing across genome
    selected_indices <- round(seq(1, n_loci, length.out = n))
    # Remove duplicates that can occur with rounding
    selected_indices <- unique(selected_indices)
    chip_vector[selected_indices] <- TRUE

  } else if (method == "chromosome_even") {
    # Proportional distribution across chromosomes
    chr_counts <- table(genome$chr)
    n_per_chr <- round(n * chr_counts / sum(chr_counts))

    # Ensure we don't exceed n due to rounding
    while (sum(n_per_chr) > n) {
      # Find chromosome with most SNPs and reduce by 1
      max_chr <- which.max(n_per_chr)
      n_per_chr[max_chr] <- n_per_chr[max_chr] - 1
    }

    # Distribute remaining SNPs if sum is less than n
    while (sum(n_per_chr) < n) {
      # Find chromosome with fewest SNPs (but still has room)
      chr_with_space <- which(n_per_chr < chr_counts)
      if (length(chr_with_space) > 0) {
        min_chr <- chr_with_space[which.min(n_per_chr[chr_with_space])]
        n_per_chr[min_chr] <- n_per_chr[min_chr] + 1
      } else {
        break
      }
    }

    for (chr in unique(genome$chr)) {
      chr_loci <- which(genome$chr == chr)
      n_chr_snp <- n_per_chr[as.character(chr)]

      if (n_chr_snp > 0 && length(chr_loci) > 0) {
        if (n_chr_snp >= length(chr_loci)) {
          # Select all loci on this chromosome
          chip_vector[chr_loci] <- TRUE
        } else {
          # Select evenly spaced loci
          selected <- round(seq(1, length(chr_loci), length.out = n_chr_snp))
          selected <- unique(selected)  # Remove duplicates from rounding
          chip_vector[chr_loci[selected]] <- TRUE
        }
      }
    }

  } else {
    stop(
      "Invalid method: '", method, "'. ",
      "Must be 'random', 'even', or 'chromosome_even'",
      call. = FALSE
    )
  }

  return(chip_vector)
}


#' Select loci by locus IDs
#'
#' @param genome Data frame of genome_meta table
#' @param locus_ids Integer vector of locus IDs to select
#' @return Logical vector indicating selected loci
#' @keywords internal
select_by_locus_ids <- function(genome, locus_ids) {

  valid_ids <- genome$locus_id

  # Check for invalid IDs
  if (!all(locus_ids %in% valid_ids)) {
    invalid <- setdiff(locus_ids, valid_ids)
    stop(
      "Invalid locus_ids: ", paste(invalid, collapse = ", "),
      call. = FALSE
    )
  }

  return(genome$locus_id %in% locus_ids)
}


#' Select loci by locus names
#'
#' @param genome Data frame of genome_meta table
#' @param locus_names Character vector of locus names to select
#' @return Logical vector indicating selected loci
#' @keywords internal
select_by_locus_names <- function(genome, locus_names) {

  valid_names <- genome$locus_name

  # Check for invalid names
  if (!all(locus_names %in% valid_names)) {
    invalid <- setdiff(locus_names, valid_names)
    stop(
      "Invalid locus_names: ", paste(head(invalid, 10), collapse = ", "),
      if (length(invalid) > 10) " ..." else "",
      call. = FALSE
    )
  }

  return(genome$locus_name %in% locus_names)
}
