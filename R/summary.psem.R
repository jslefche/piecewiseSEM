#' Evaluate a list of structural equations with grouping
#'
#' @param modelList a list of structural equations

summary.psem <- function(modelList, groups = NULL,
                         direction = NULL, conserve = FALSE, conditional = FALSE,
                         add.claims = NULL,
                         intercepts = FALSE, standardize = TRUE,
                         .progressBar = TRUE) {

  if(!is.null(groups)) {

    data <- modelList$data

    lvls <- unique(data[, colnames(data) %in% groups])

    modelList2 <- lapply(lvls, function(i) update(modelList, data = data[data[, colnames(data) %in% groups] == i, ]))

    names(modelList2) <- lvls

    lapply(modelList2, function(i) summary.psem2(i, direction, conserve, conditional, add.claims, intercepts, standardize, .progressBar))

  } else summary.psem2(modelList, direction, conserve, conditional, add.claims, intercepts, standardize, .progressBar)

}

#' Evaluate a list of structural equations
summary.psem2 <- function(modelList,
                         direction = NULL, conserve = FALSE, conditional = FALSE,
                         add.claims = NULL,
                         intercepts = FALSE, standardize = TRUE,
                         .progressBar = TRUE) {

  name <- deparse(substitute(modelList))

  call <- paste(listFormula(modelList), collapse = "\n  ")

  dTable <- dSep(modelList, direction, conserve, conditional, .progressBar)

  Cstat <- fisherC(dTable, add.claims, direction, conserve, conditional, .progressBar)

  IC <- infCrit(modelList, Cstat, add.claims, direction, conserve, conditional, .progressBar)

  coefficients <- coefs(modelList, intercepts, standardize)

  R2 <- rsquared(modelList)

  R2[, which(sapply(R2, is.numeric))] <- round(R2[, which(sapply(R2, is.numeric))], 2)

  if(length(dTable) > 0)

    dTable[, which(sapply(dTable, is.numeric))] <- round(dTable[, which(sapply(dTable, is.numeric))], 4)

  l <- list(name = name, call = call, dTable = dTable, Cstat = Cstat, IC = IC, coefficients = coefficients, R2 = R2)

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

  cat("\n---\nTests of directed separation:\n\n", captureTable(x$dTable))

  cat("\nGlobal goodness-of-fit:\n\n  Fisher's C =", as.character(x$Cstat[1]),
      "with P-value =", as.character(x$Cstat[3]),
      "and on", as.character(x$Cstat[2]), "degrees of freedom")

  cat("\n\n---\nCoefficients:\n\n", captureTable(x$coefficients))

  cat("\n  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘’ 1")

  cat("\n\nIndividual R-squared:\n\n", captureTable(x$R2[, c(1, 4:ncol(x$R2))]))

  invisible(x)

}

captureTable <- function(g) {

  g1 <- capture.output(print(g, row.names = FALSE))

  g1 <- paste0(g1, "\n")

  return(g1)

}
