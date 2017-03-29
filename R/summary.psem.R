#' Evaluate a list of structural equations
#'
#' @param modelList a list of structural equations

summary.psem <- function(modelList, direction = NULL, conserve = FALSE, conditional = FALSE,
                         intercepts = FALSE, standardize = TRUE,
                         .progressBar = TRUE) {

  if(class(modelList) == "psem") data <- modelList$data

  if(is.null(data)) data <- getData.(modelList[[1]])

  name <- deparse(substitute(modelList))

  call <- paste(listFormula(modelList), collapse = "\n  ")

  dTable <- dSep(modelList, direction, conserve, conditional, .progressBar)

  C <- fisherC(dTable)

  IC <- InfCrit(modelList, C)

  coefficients <- coefs(modelList, data, intercepts, standardize)

  # R2 <- rsquared(modelList)

  if(length(dTable) > 0)

    dTable[, which(sapply(dTable, is.numeric))] <- round(dTable[, which(sapply(dTable, is.numeric))], 4)

  l <- list(name = name, call = call, dTable = dTable, C = C, IC = IC, coefficients = coefficients) #, R2 = R2)

  class(l) <- "summary.psem"

  l

}

print.summary.psem <- function(x) {

  cat("\nStructural Equation Model of", as.character(x$name), "\n")

  cat("\nCall:\n ", x$call)

  cat("\n")

  cat("\n    AIC      BIC")
  cat("\n", as.character(sprintf("%.3f", x$IC[1])), " ", as.character(x$IC[3]))

  cat("\n")

  cat("\nTests of directed separation:\n\n", captureTable(x$dTable))

  cat("\nCoefficients:\n\n", captureTable(x$coefficients))

  cat("  ---\nSignif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘’ 1")

  cat("\n\n")

  cat("Goodness-of-fit:\n\n  Global model: Fisher's C =", as.character(x$C[1]),
      "with P-value =", as.character(x$C[3]),
      "and on", as.character(x$C[2]), "degrees of freedom")

  # cat("\n  Individual models: R-squared =\n  ", captureTable(x$R2))

  invisible(x)

}

captureTable <- function(g) {

  g1 <- capture.output(print(g, row.names = FALSE))

  g1 <- paste0(g1, "\n")

  return(g1)

}
