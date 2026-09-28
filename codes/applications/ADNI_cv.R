# Fixed design: no outcome or fitted performance is used to choose the split.
ADNI_CV_FOLDS <- 5L
ADNI_CV_SEED <- 20260914L
adni_make_folds <- function(dat) {
  stopifnot(nrow(dat) == 495L, !anyNA(dat$RID), !anyDuplicated(dat$RID),
            !anyNA(dat$R), all(dat$R %in% 0:1), sum(dat$R == 1L) == 133L)
  old_kind <- RNGkind()
  old_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({ do.call(RNGkind, as.list(old_kind))
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(ADNI_CV_SEED)
  z <- dat[order(dat$RID), c("RID", "R")]
  z$fold <- NA_integer_; offset <- 0L
  # Continuing the cycle across strata gives 99 people in every fold and
  # either 26 or 27 observed-CSF participants in each held-out group.
  for (r in 0:1) {
    idx <- which(z$R == r)
    idx <- idx[sample.int(length(idx))]
    z$fold[idx] <- (offset + seq_along(idx) - 1L) %% ADNI_CV_FOLDS + 1L
    offset <- offset + length(idx)
  }
  stopifnot(all(table(z$fold) == 99L), all(table(z$fold[z$R == 1L]) %in% 26:27))
  rownames(z) <- NULL
  z
}
adni_read_cohort <- function(data_dir) {
  path <- file.path(data_dir, "toledo_adni_cohort_n495.csv")
  if (unname(tools::md5sum(path)) != "d10b980253304620efea9f2a013166a9")
    stop("The cleaned cohort differs from the frozen ADNI cohort: ", path)
  dat <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("RID", "R", "y_centiloid", "u_csf_abeta42", "x_age_years",
                "x_female", "x_apoe4_dose", "c_ab42_ab40_3", "c_gfap_3", "diagnosis")
  stopifnot(all(required %in% names(dat)))
  dat
}
