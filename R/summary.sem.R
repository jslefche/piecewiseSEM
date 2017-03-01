#' Evaluate a list of structural equations

#' @param modelList a list of structural equations

summary.sem <- function(modelList, conditional = FALSE, .progressBar = TRUE) {

  name <- deparse(substitute(modelList))

  call <- paste(listFormula(modelList), collapse = "\n  ")

  dTable <- dSep(modelList, conditional, .progressBar)

  C <- fisherC(dTable)

  IC <- InfCrit(modelList, C)

  l <- list(name = name, call = call, dTable = dTable, C = C, IC = IC)

  class(l) <- "summary.sem"

  l

}

print.summary.sem <- function(x) {

  cat("\nStructural Equation Model of", as.character(x$name), "\n")

  cat("\nCall:\n ", x$call)

  cat("\n")

  cat("\n    AIC      BIC")
  cat("\n", as.character(sprintf("%.3f", x$IC[1])), " ", as.character(x$IC[3]))

  cat("\n")

  cat("\nTests of directed separation:\n\n", captureTable(print.data.frame(x$dTable, row.names = FALSE)))

  cat("---\nSignif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘’ 1")

  cat("\n\n\n")

  cat("Goodness-of-fit:\n\n Fisher's C =", as.character(x$C[1]),
      "with P-value =", as.character(x$C[3]),
      "and on", as.character(x$C[2]), "degrees of freedom")

  cat("\n\n")

  invisible(x)

}

captureTable <- function(x) {

  x <- capture.output(x)

  x <- paste0(x, "\n")

  return(x)

}
