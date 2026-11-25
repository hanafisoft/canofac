#' @export
print.canofac <- function(x, digits = 3, ...) {
  s <- x$settings
  cat(.tr("title_fit"), "\n", sep = "")
  cat(sprintf(.tr("settings_fit"), s$q, s$axes, s$extraction, s$R, nrow(x$d)))
  if (!is.null(s$weighting) && s$weighting != "none")
    cat(sprintf(.tr("weighting_line"), toupper(s$weighting)), "\n", sep = "")
  cat(sprintf(.tr("criterion"), format(x$criterion, digits = 6)))
  cat(.tr("pvs_head"), "\n", sep = "")
  print(round(x$d, digits))
  invisible(x)
}

#' @export
summary.canofac <- function(object, ...) {
  out <- list(d = object$d, saliences = object$saliences,
              delta = object$delta, criterion = object$criterion,
              settings = object$settings)
  print(object, ...)
  cat("\n", .tr("sal_head"), "\n", sep = "")
  print(round(object$saliences, 3))
  cat("\n", .tr("delta_head"), "\n", sep = "")
  print(round(object$delta, 3))
  invisible(out)
}

#' @export
print.cf_test <- function(x, digits = 1, ...) {
  s <- x$settings
  fixed <- s$statistic %in% c("lobo", "fixed")
  cat(.tr("title_test"), "\n", sep = "")
  cat(sprintf(.tr("settings_test"), s$q, s$axes, s$extraction, s$scheme, s$B, s$alpha))
  cat(.tr(if (fixed) "stat_fixed" else "stat_refit"), "\n", sep = "")
  if (fixed) cat(.tr(if (identical(s$subspace, "inertia")) "sub_inertia" else "sub_equal"), "\n", sep = "")
  if (!is.null(s$weighting) && s$weighting != "none")
    cat(sprintf(.tr("weighting_line"), toupper(s$weighting)), "\n", sep = "")
  if (s$scheme %in% c("groups", "within"))
    cat(sprintf(.tr("groups_line"), s$n_groups, s$scheme), "\n", sep = "")
  pv <- if (is.null(s$pvalue)) "permutation" else s$pvalue
  if (pv == "permutation") {
    cat(sprintf(.tr("pvalue_perm"), s$B), "\n", sep = "")
  } else if (pv == "gamma") {
    cat(.tr("pvalue_gamma"), "\n", sep = "")
  } else {
    cat(sprintf(.tr("pvalue_hybrid"), s$B, s$hybrid_c, s$hybrid_c,
                sum(x$pmethod == "perm"), length(x$pmethod)), "\n", sep = "")
  }
  cat("\n")
  cat(.tr("psub_head"), "\n", sep = "")
  print(round(100 * x$p_subspace, digits))
  who <- names(x$p_subspace)[x$linked_subspace]
  cat(sprintf(.tr("subspace_line"), .link_label(x$link_subspace),
              if (length(who)) paste0(" (", paste(who, collapse = ", "), ")") else ""), "\n", sep = "")
  cat("\n", .tr(if (fixed) "pval_head_fixed" else "pval_head"), "\n", sep = "")
  print(round(100 * x$p, digits))
  cat("\n", .tr(if (fixed) "struct_head_fixed" else "struct_head"), "\n", sep = "")
  lab <- .link_label(x$link)
  w <- max(nchar(lab))
  for (h in seq_along(x$link)) {
    who <- rownames(x$linked)[x$linked[, h]]
    cat(sprintf("  %s: %-*s %s\n", names(x$link)[h], w, lab[h],
                if (length(who)) paste0("(", paste(who, collapse = ", "), ")") else ""))
  }
  yn <- function(v) ifelse(v, .tr("yes"), .tr("no"))
  cat("\n", .tr(if (fixed) "linked_line" else "retained_line"),
      paste(names(x$retained), yn(x$retained), sep = "=", collapse = ", "), "\n", sep = "")
  if (!is.null(x$p_compare)) {
    cat("\n", .tr("compare_head"), "\n", sep = "")
    print(round(100 * x$p_compare, digits))
    cat(.tr("compare_note"), "\n", sep = "")
  }
  invisible(x)
}
