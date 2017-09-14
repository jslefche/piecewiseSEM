#' Evaluate a list of structural equations with grouping
#'
#' @param object a list of structural equations

summary.psem <- function(object, groups = NULL,
                         direction = NULL, conserve = FALSE, conditional = FALSE,
                         add.claims = NULL,
                         intercepts = FALSE, standardize = TRUE,
                         .progressBar = TRUE) {

  if(!is.null(groups)) {

    data <- object$data

    lvls <- unique(data[, colnames(data) %in% groups])

    object2 <- lapply(lvls, function(i) update(object, data = data[data[, colnames(data) %in% groups] == i, ]))

    names(object2) <- lvls

    lapply(object2, function(i) summary2.psem(i, direction, conserve, conditional, add.claims, intercepts, standardize, .progressBar))

  } else summary2.psem(object, direction, conserve, conditional, add.claims, intercepts, standardize, .progressBar)

}

#' Evaluate a list of structural equations
summary2.psem <- function(object,
                         direction = NULL, conserve = FALSE, conditional = FALSE,
                         add.claims = NULL,
                         intercepts = FALSE, standardize = TRUE,
                         .progressBar = TRUE) {

  name <- deparse(substitute(object))

  call <- paste(listFormula(object), collapse = "\n  ")

  dTable <- dSep(object, direction, conserve, conditional, .progressBar)

  Cstat <- fisherC(dTable, add.claims, direction, conserve, conditional, .progressBar)

  IC <- infCrit(object, Cstat, add.claims, direction, conserve, conditional, .progressBar)

  coefficients <- coefs(object, intercepts, standardize)

  R2 <- rsquared(object)

  R2[, which(sapply(R2, is.numeric))] <- round(R2[, which(sapply(R2, is.numeric))], 2)

  if(length(dTable) > 0)

    dTable[, which(sapply(dTable, is.numeric))] <- round(dTable[, which(sapply(dTable, is.numeric))], 4)

  l <- list(name = name, call = call, dTable = dTable, Cstat = Cstat, IC = IC, coefficients = coefficients, R2 = R2)

  class(l) <- "summary.psem"

  l

}

print.summary.psem <- function(object, ...) {

  cat("\nStructural Equation Model of", as.character(object$name), "\n")

  cat("\nCall:\n ", object$call)

  cat("\n")

  cat("\n    AIC      BIC")
  cat("\n", as.character(sprintf("%.3f", object$IC[1])), " ", as.character(object$IC[3]))

  cat("\n")

  cat("\n---\nTests of directed separation:\n\n", captureTable(object$dTable))

  cat("\nGlobal goodness-of-fit:\n\n  Fisher's C =", as.character(object$Cstat[1]),
      "with P-value =", as.character(object$Cstat[3]),
      "and on", as.character(object$Cstat[2]), "degrees of freedom")

  cat("\n\n---\nCoefficients:\n\n", captureTable(object$coefficients))

  cat("\n  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05")

  cat("\n\nIndividual R-squared:\n\n", captureTable(object$R2[, c(1, 4:ncol(object$R2))]))

  invisible(object)

}

captureTable <- function(g) {

  g1 <- capture.output(print(g, row.names = FALSE))

  g1 <- paste0(g1, "\n")

  return(g1)

}
