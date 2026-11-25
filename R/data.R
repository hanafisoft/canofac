#' Friday (1987) macro-invertebrate data
#'
#' Loads the Friday (1987) data distributed with the \pkg{ade4} package (16
#' ponds, 91 species in 10 faunistic groups) and returns the 10 blocks,
#' centred by column (species are not scaled) and at unit total inertia, as in
#' Chessel and Hanafi (1996).
#'
#' @param preprocess logical, apply \code{\link{cf_preprocess}} (centring and
#'   unit total inertia).
#' @return A named list of 10 matrices (16 rows).
#' @references Friday L.E. (1987). The diversity of macroinvertebrate and
#'   macrophyte communities in ponds. Freshwater Biology, 18, 87--104.
#'
#'   Chessel D., Hanafi M. (1996). Analyses de la co-inertie de K nuages de
#'   points. Revue de Statistique Appliquée, 44(2), 35--60.
#' @export
cf_friday <- function(preprocess = TRUE) {
  if (!requireNamespace("ade4", quietly = TRUE))
    stop("package 'ade4' is required: install.packages('ade4')")
  env <- new.env()
  utils::data("friday87", package = "ade4", envir = env)
  f <- env$friday87
  fau <- as.matrix(f$fau)
  if (nrow(fau) != 16 && ncol(fau) == 16) fau <- t(fau)
  blo <- f$fau.blo
  grp <- rep(seq_along(blo), blo)
  blocks <- lapply(seq_along(blo), function(g) fau[, grp == g, drop = FALSE])
  names(blocks) <- as.character(f$tab.names)
  if (preprocess) blocks <- cf_preprocess(blocks)
  blocks
}

#' breast.TCGA training data
#'
#' Loads the training data of \code{breast.TCGA} from \pkg{mixOmics} (150
#' tumours; mRNA, miRNA and proteins). The data were normalised and filtered
#' by the authors of the package; here each variable is centred and scaled and
#' each block is brought to unit total inertia.
#'
#' @param preprocess logical, centre, scale and bring each block to unit total
#'   inertia.
#' @return A list with \code{blocks} (named list of 3 matrices) and
#'   \code{subtype} (factor).
#' @references Rohart F., Gautier B., Singh A., Lê Cao K.-A. (2017). mixOmics:
#'   an R package for 'omics feature selection and multiple data integration.
#'   PLoS Computational Biology, 13(11), e1005752.
#' @export
cf_tcga <- function(preprocess = TRUE) {
  if (!requireNamespace("mixOmics", quietly = TRUE))
    stop("package 'mixOmics' is required (Bioconductor)")
  env <- new.env()
  utils::data("breast.TCGA", package = "mixOmics", envir = env)
  tr <- env$breast.TCGA$data.train
  blocks <- list(mRNA = tr$mrna, miRNA = tr$mirna, protein = tr$protein)
  if (preprocess) blocks <- cf_preprocess(blocks, scale = TRUE)
  list(blocks = blocks, subtype = tr$subtype)
}

#' Potato data (articles A and BB)
#'
#' Potato data of Thygesen et al. (2001), distributed with the R package
#' multiblock: 26 samples, 9 blocks. Two sets of blocks are predefined:
#' \code{"BB"} (article on the LOBO test) keeps Chemical (14 variables),
#' Compression (12) and Sensory (9); \code{"A"} (article on pseudo-singular
#' values) keeps the three blocks measured on the raw tubers, NIRraw (1050
#' wavelengths), FIDraw (500 points) and CPMGraw (410 points). In the
#' multiblock version, each NIR spectrum is SNV-transformed and the variables
#' of Chemical and Compression are already scaled; preprocessing: centring and
#' unit total inertia.
#'
#' @param file optional CSV file in the exchange format of
#'   \code{\link{cf_read_blocks}}; if \code{NULL}, the data are read from the
#'   package multiblock.
#' @param blocks blocks to keep: a character vector of block names, or
#'   \code{"BB"} or \code{"A"} for the predefined sets.
#' @param preprocess logical: centre and bring each block to unit total inertia.
#' @return A named list of matrices.
#' @references Thygesen L. G., Thybo A. K., Engelsen S. B. (2001). Prediction
#'   of sensory texture quality of boiled potatoes from low-field 1H NMR of raw
#'   potatoes. The role of chemical constituents. LWT - Food Science and
#'   Technology, 34, 469-477.
#'
#'   Liland K. H., Naes T., Indahl U. G. (2016). ROSA - a fast extension of
#'   partial least squares regression for multiblock data analysis. Journal of
#'   Chemometrics, 30, 651-662.
#' @examples
#' \dontrun{
#' blocks <- cf_potato()            # article BB
#' cf_test(blocks, R = 2, q = 2, B = 999, seed = 1)
#' raw <- cf_potato(blocks = "A")   # article A
#' canofac(raw, R = 2, q = 4)
#' }
#' @export
cf_potato <- function(file = NULL, blocks = "BB", preprocess = TRUE) {
  sets <- list(BB = c("Chemical", "Compression", "Sensory"),
               A = c("NIRraw", "FIDraw", "CPMGraw"))
  if (length(blocks) == 1 && blocks %in% names(sets)) blocks <- sets[[blocks]]
  if (!is.null(file)) {
    bl <- cf_read_blocks(file, blocks)$blocks
  } else {
    if (!requireNamespace("multiblock", quietly = TRUE))
      stop("package 'multiblock' is required, or give a file")
    env <- new.env()
    utils::data("potato", package = "multiblock", envir = env)
    pot <- env$potato
    miss <- setdiff(blocks, names(pot))
    if (length(miss)) stop("blocks not found: ", paste(miss, collapse = ", "))
    bl <- lapply(blocks, function(b) {
      X <- as.matrix(pot[[b]])
      storage.mode(X) <- "double"
      X
    })
    names(bl) <- blocks
  }
  if (preprocess) bl <- cf_preprocess(bl)
  bl
}

#' Read and write multiblock data (exchange format shared with Python)
#'
#' One CSV file: first column the row names, then the variables named
#' \code{<block>.x.<variable>} (layout of the potato data of the package
#' multiblock). Columns without \code{.x.} are extra variables (for example a
#' group). This format is read and written by the R and Python versions of
#' canofac, for example to export the datasets of ade4, mixOmics or multiblock
#' for Python: \code{cf_write_blocks(cf_friday(preprocess = FALSE), "friday.csv")}.
#'
#' @param file CSV file.
#' @param blocks for \code{cf_read_blocks}, blocks to keep (default all, in file
#'   order); for \code{cf_write_blocks}, a named list of matrices or data frames.
#' @param extra optional data frame of extra columns to write.
#' @return \code{cf_read_blocks}: a list with \code{blocks} (named list of
#'   matrices) and \code{extra} (data frame of the other columns).
#' @examples
#' f <- tempfile(fileext = ".csv")
#' bl <- list(A = matrix(rnorm(20), 5), B = matrix(rnorm(15), 5))
#' cf_write_blocks(bl, f)
#' str(cf_read_blocks(f)$blocks)
#' @export
cf_read_blocks <- function(file, blocks = NULL) {
  d <- utils::read.csv(file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  cn <- colnames(d)
  isb <- grepl(".x.", cn, fixed = TRUE)
  bname <- ifelse(isb, sub("\\.x\\..*$", "", cn), NA)
  vname <- ifelse(isb, sub("^.*?\\.x\\.", "", cn, perl = TRUE), NA)
  all_b <- unique(bname[isb])
  keep <- if (is.null(blocks)) all_b else blocks
  miss <- setdiff(keep, all_b)
  if (length(miss)) stop("blocks not found: ", paste(miss, collapse = ", "))
  out <- lapply(keep, function(b) {
    j <- which(bname == b)
    X <- as.matrix(d[, j, drop = FALSE])
    colnames(X) <- vname[j]
    storage.mode(X) <- "double"
    X
  })
  names(out) <- keep
  list(blocks = out, extra = d[, !isb, drop = FALSE])
}

#' @rdname cf_read_blocks
#' @export
cf_write_blocks <- function(blocks, file, extra = NULL) {
  parts <- lapply(names(blocks), function(b) {
    X <- as.data.frame(blocks[[b]], check.names = FALSE)
    if (is.null(colnames(blocks[[b]]))) colnames(X) <- seq_len(ncol(X)) - 1
    colnames(X) <- paste0(b, ".x.", colnames(X))
    X
  })
  out <- do.call(cbind, parts)
  rn <- rownames(blocks[[1]])
  rownames(out) <- if (is.null(rn)) seq_len(nrow(out)) - 1 else rn
  if (!is.null(extra)) out <- cbind(out, extra)
  utils::write.csv(out, file)
  invisible(file)
}
