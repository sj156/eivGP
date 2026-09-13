#!/usr/bin/env Rscript

## Build the private assignment files for the 20-validation ADNI protocol.
##
## Validations 1--9 are copied exactly from the frozen balanced 3 x 3 file.
## Validations 10--20 are separately randomized 1/3 test holdouts within
## observed-CSF status only.  No outcome, observed CSF value, diagnosis, or
## fitted-model performance is used to choose the 11 new assignments.

options(stringsAsFactors = FALSE)

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Run this file with Rscript.")
script_file <- normalizePath(sub("^--file=", "", file_arg))

script_dir <- dirname(script_file)
project_dir <- dirname(script_dir)

source(file.path(project_dir, "..", "adni_paths.R"), local = TRUE)
data_dir <- resolve_adni_data_dir(project_dir)

output_dir <- Sys.getenv(
  "ADNI_SPLIT_DIR",
  unset = file.path(data_dir, "validation-splits")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#data_dir <- resolve_adni_data_dir(project_dir)
#output_dir <- Sys.getenv(
#  "ADNI_SPLIT_DIR", unset = file.path(data_dir, "validation-splits")
#)
#dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cohort_file <- file.path(data_dir, "toledo_adni_cohort_n495.csv")
legacy_file <- file.path(data_dir, "toledo_adni_balanced_repeated_3fold.csv")
plan_file <- file.path(project_dir, "validation_split_plan.csv")
expected_cohort_md5 <- "d10b980253304620efea9f2a013166a9"
expected_legacy_md5 <- "895a5759c7baf22e1771516c0aa84d70"

stopifnot(
  unname(tools::md5sum(cohort_file)) == expected_cohort_md5,
  unname(tools::md5sum(legacy_file)) == expected_legacy_md5,
  file.exists(plan_file)
)

atomic_write_csv <- function(x, path) {
  temporary <- paste0(path, ".tmp")
  utils::write.csv(x, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) stop("Could not atomically write ", path)
}
bin_rank <- function(x, k) {
  pmin(k, pmax(1L, ceiling(rank(x, ties.method = "average") / length(x) * k)))
}
standardized_difference <- function(test, full) {
  scale <- stats::sd(full)
  if (!is.finite(scale) || scale <= 0) return(NA_real_)
  (mean(test) - mean(full)) / scale
}

cohort <- utils::read.csv(cohort_file, check.names = FALSE)
legacy <- utils::read.csv(legacy_file, check.names = FALSE)
plan <- utils::read.csv(plan_file, check.names = FALSE)
stopifnot(
  nrow(cohort) == 495L, !anyDuplicated(cohort$RID), sum(cohort$R == 1L) == 133L,
  nrow(legacy) == 1485L, !anyDuplicated(legacy[, c("repeat_id", "RID")]),
  identical(plan$validation_id, 1:20),
  all(plan$test_n == 165L), all(plan$test_R1_n %in% 44:45)
)
cohort <- cohort[order(cohort$RID), , drop = FALSE]
r1_rows <- which(cohort$R == 1L)
r0_rows <- which(cohort$R == 0L)
cell <- (cohort$c_ab42_ab40_3 - 1L) * 3L + cohort$c_gfap_3

## A support condition prevents a random split from making an observed-CSF
## proxy cell unusable in training.  This is a prespecified design check, not
## a predictive-performance gate.  It is the only rejection rule.
draw_random_holdout <- function(seed, test_R1_n, max_attempts = 10000L) {
  set.seed(seed)
  for (attempt in seq_len(max_attempts)) {
    test_rows <- c(
      sample(r1_rows, test_R1_n, replace = FALSE),
      sample(r0_rows, 165L - test_R1_n, replace = FALSE)
    )
    train_r1 <- setdiff(r1_rows, test_rows)
    support <- table(factor(cell[train_r1], levels = 1:9))
    if (min(support) >= 2L) {
      return(list(test_rows = sort(test_rows), rejection_attempt = attempt))
    }
  }
  stop("No support-valid split found for seed ", seed)
}

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
assignments <- vector("list", 20L)
audit <- vector("list", 20L)
file_manifest <- vector("list", 20L)
test_signatures <- character(20L)

for (i in seq_len(nrow(plan))) {
  row <- plan[i, , drop = FALSE]
  validation_id <- row$validation_id
  if (row$source == "legacy_balanced_3x3") {
    old <- legacy[
      legacy$repeat_id == row$legacy_repeat_id &
        legacy$fold == row$legacy_fold_id, , drop = FALSE
    ]
    test_rows <- match(old$RID, cohort$RID)
    if (anyNA(test_rows)) stop("Legacy split contains an unknown RID.")
    rejection_attempt <- NA_integer_
    generation_seed <- row$assignment_seed
  } else if (row$source == "random_holdout_v1") {
    answer <- draw_random_holdout(row$sampling_seed, row$test_R1_n)
    test_rows <- answer$test_rows
    rejection_attempt <- answer$rejection_attempt
    generation_seed <- row$sampling_seed
  } else stop("Unknown source in validation_split_plan.csv: ", row$source)

  is_test <- seq_len(nrow(cohort)) %in% test_rows
  if (sum(is_test) != row$test_n || sum(cohort$R[is_test] == 1L) != row$test_R1_n) {
    stop("Validation ", validation_id, " failed its frozen size/R1 check.")
  }
  test_signatures[i] <- paste(sort(cohort$RID[is_test]), collapse = ",")
  split <- data.frame(
    validation_id = validation_id,
    RID = cohort$RID,
    partition = ifelse(is_test, "test", "train"),
    source = row$source,
    generation_seed = generation_seed,
    rejection_attempt = rejection_attempt,
    legacy_repeat_id = row$legacy_repeat_id,
    legacy_fold_id = row$legacy_fold_id
  )
  split_file <- file.path(output_dir, sprintf("validation_%02d.csv", validation_id))
  atomic_write_csv(split, split_file)
  assignments[[i]] <- split

  test <- cohort[is_test, , drop = FALSE]
  train_r1 <- cohort[!is_test & cohort$R == 1L, , drop = FALSE]
  test_r1 <- test[test$R == 1L, , drop = FALSE]
  support <- table(factor(
    (train_r1$c_ab42_ab40_3 - 1L) * 3L + train_r1$c_gfap_3,
    levels = 1:9
  ))
  audit[[i]] <- data.frame(
    validation_id = validation_id, source = row$source,
    test_n = nrow(test), train_n = nrow(cohort) - nrow(test),
    test_R1 = nrow(test_r1), train_R1 = nrow(train_r1),
    age_SMD = standardized_difference(test$x_age_years, cohort$x_age_years),
    Y_SMD = standardized_difference(test$y_centiloid, cohort$y_centiloid),
    U_SMD_R1 = standardized_difference(
      test_r1$u_csf_abeta42, cohort$u_csf_abeta42[cohort$R == 1L]
    ),
    min_train_R1_Ccell = min(support),
    generation_seed = generation_seed, rejection_attempt = rejection_attempt
  )
  file_manifest[[i]] <- data.frame(
    validation_id = validation_id, file = basename(split_file),
    md5 = unname(tools::md5sum(split_file)), source = row$source,
    generation_seed = generation_seed, rejection_attempt = rejection_attempt
  )
}

if (anyDuplicated(test_signatures)) stop("Two validation test sets are identical.")
all_assignments <- do.call(rbind, assignments)
all_audit <- do.call(rbind, audit)
manifest <- do.call(rbind, file_manifest)
stopifnot(
  nrow(all_assignments) == 20L * nrow(cohort),
  !anyDuplicated(all_assignments[, c("validation_id", "RID")]),
  all(table(all_assignments$validation_id, all_assignments$partition)[, "test"] == 165L),
  all(all_audit$min_train_R1_Ccell >= 2L)
)

atomic_write_csv(
  all_assignments,
  file.path(output_dir, "all_validation_assignments.csv")
)
atomic_write_csv(all_audit, file.path(output_dir, "validation_balance_audit.csv"))
atomic_write_csv(manifest, file.path(output_dir, "validation_file_manifest.csv"))

test_sets <- lapply(assignments, function(x) x$RID[x$partition == "test"])
overlap <- do.call(rbind, lapply(1:20, function(i) {
  do.call(rbind, lapply(1:20, function(j) data.frame(
    validation_i = i, validation_j = j,
    n_test_intersection = length(intersect(test_sets[[i]], test_sets[[j]])),
    Jaccard = length(intersect(test_sets[[i]], test_sets[[j]])) /
      length(union(test_sets[[i]], test_sets[[j]]))
  )))
}))
atomic_write_csv(overlap, file.path(output_dir, "validation_test_overlap.csv"))

writeLines(c(
  "Private ADNI 20-validation assignments",
  "",
  "validation_01.csv through validation_20.csv contain RID and train/test membership.",
  "Validations 1--9 exactly copy the frozen balanced 3 x 3 assignments.",
  "Validations 10--20 are R-stratified random 1/3 holdouts with no outcome-based selection.",
  "The only rejection condition is at least two observed-CSF training participants in every 3 x 3 proxy cell.",
  "These files contain participant identifiers. Transfer them only through approved private storage; never commit them to the public repository.",
  "The 20 validations overlap; across-validation SD is a descriptive stability summary."
), file.path(output_dir, "README_PRIVATE.txt"))

cat("Wrote 20 private validation assignments to:", output_dir, "\n")
cat("Combined assignment MD5:", unname(tools::md5sum(
  file.path(output_dir, "all_validation_assignments.csv")
)), "\n")
print(all_audit, row.names = FALSE)
