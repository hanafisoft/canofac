#' Language of the printed and plotted results
#'
#' Sets or returns the language used for the labels of the results: printed
#' summaries (\code{print}, \code{summary}), plot titles, axis labels and
#' legends, and the outcomes of the test (all blocks, several blocks, one
#' block, no block linked). The values stored in the objects, for example
#' \code{link} in a \code{cf_test} object, are not translated, so that
#' code that uses them does not depend on the language.
#'
#' The language is stored in the option \code{canofac.lang}; it can also be set
#' with \code{options(canofac.lang = "fr")}, for instance in a
#' \code{.Rprofile} file. The labels are read from the file
#' \code{inst/shared/labels.csv}, shared with the Python version of canofac.
#' The help pages exist in English (\code{?cf_test}) and in French
#' (\code{cf_help("cf_test", "fr")}).
#'
#' @param lang \code{"en"} (English, default) or \code{"fr"} (French). If
#'   missing, the current language is returned.
#' @return The current language, invisibly when it is set.
#' @examples
#' cf_lang()            # current language
#' old <- cf_lang("fr") # French labels
#' cf_lang(old)         # back to the previous language
#' @export
cf_lang <- function(lang) {
  cur <- .lang()
  if (missing(lang)) return(cur)
  lang <- match.arg(lang, c("en", "fr"))
  options(canofac.lang = lang)
  invisible(cur)
}

.lang <- function() {
  l <- getOption("canofac.lang", "en")
  if (!l %in% c("en", "fr")) "en" else l
}

# Dictionary of labels: read once from inst/shared/labels.csv, the file also
# used by the Python version, so that both print and plot the same words.
.labels_env <- new.env(parent = emptyenv())

.dict_table <- function() {
  if (is.null(.labels_env$d)) {
    f <- system.file("shared", "labels.csv", package = "canofac")
    x <- utils::read.csv(f, stringsAsFactors = FALSE, encoding = "UTF-8",
                         fileEncoding = "UTF-8", na.strings = character(0))
    x$en <- gsub("\\n", "\n", x$en, fixed = TRUE)
    x$fr <- gsub("\\n", "\n", x$fr, fixed = TRUE)
    .labels_env$d <- x
  }
  .labels_env$d
}

# Label for a key in the current language.
.tr <- function(key) {
  d <- .dict_table()
  i <- match(key, d$key)
  if (is.na(i)) return(key)
  if (.lang() == "fr") d$fr[i] else d$en[i]
}

# Translated outcomes of the test ("all", "several", "one", "none").
.link_label <- function(link) {
  vapply(as.character(link), .tr, character(1), USE.NAMES = FALSE)
}
