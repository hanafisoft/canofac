#' Help in English or French
#'
#' In English, the standard help page (\code{?cf_test}). In French, the pages
#' of \code{inst/shared/aide/fr}, shared with the Python version of canofac.
#' Without \code{lang}, the language of \code{\link{cf_lang}} is used.
#'
#' @param topic function name (\code{"cf_test"}, ...), \code{"canofac-package"}
#'   for the overview, or \code{"topics"} to list the French topics.
#' @param lang \code{"en"} or \code{"fr"}.
#' @return Invisibly, the text of the French page (or \code{NULL}).
#' @examples
#' cf_help("topics")
#' cf_help("cf_test", "fr")
#' @export
cf_help <- function(topic = "canofac-package", lang = NULL) {
  if (is.null(lang)) lang <- .lang()
  dir <- system.file("shared", "aide", "fr", package = "canofac")
  topics <- sub("\\.md$", "", list.files(dir, pattern = "\\.md$"))
  if (identical(topic, "topics")) {
    cat(paste(topics, collapse = ", "), "\n")
    return(invisible(NULL))
  }
  if (lang == "fr") {
    f <- file.path(dir, paste0(topic, ".md"))
    if (!file.exists(f)) {
      cat(sprintf(.tr("help_missing"), topic, paste(topics, collapse = ", ")), "\n")
      return(invisible(NULL))
    }
    txt <- readLines(f, encoding = "UTF-8", warn = FALSE)
    cat(txt, sep = "\n")
    return(invisible(txt))
  }
  print(utils::help(topic, package = "canofac"))
  invisible(NULL)
}
