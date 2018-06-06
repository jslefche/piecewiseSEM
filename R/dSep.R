#' Tests of directed separation
#' 
#' Evaluation of conditional independence claims to be used in determining the
#' goodness-of-fit for piecewise structural equation models.
#' 
#' In cases involving non-normally distributed responses in the independence
#' claims that are modeled using generalized linear models, the significance of
#' the independence claim is not reversable (e.g., the P-value of Y ~ X is not
#' the same as X ~ Y). This is due to the transformation of the response via
#' the link function. In extreme cases, this can bias the goodness-of-fit
#' tests. \code{summary.psem} will issue a warning when this case is present
#' and provide guidance for solutions.
#' 
#' One solution is to specify the directionality of the relationship using the
#' \code{direction} argument, e.g. \code{direction = c("X <- Y")}. Another is
#' to run both tests (Y ~ X, X ~ Y) and return the most conservative (i.e.,
#' lowest) P-value, which can be toggled using the \code{conserve = TRUE}
#' argument.
#' 
#' @param modelList A list of structural equations created using \code{psem}.
#' @param direction A \code{vector} of claims defining the specific
#' directionality of independence claims; for use in special cases (see
#' Details).
#' @param conserve Whether the most conservative P-value should be returned;
#' for use in special cases (see Details). Default is FALSE.
#' @param conditioning Whether the conditioning variables should be shown in
#' the summary table. Default is FALSE.
#' @param .progressBar An optional progress bar. Default is TRUE.
#' @return Returns a \code{data.frame} of independence claims and their
#' significance values.
#' @author Jon Lefcheck <jlefcheck@@bigelow.org>
#' @seealso \code{\link{basisSet}}
#' @references Shipley, Bill. "A new inferential test for path models based on
#' directed acyclic graphs." Structural Equation Modeling 7.2 (2000): 206-218.
#' @examples
#' 
#' 
#' 
#' @export dSep
dSep <- function(modelList, direction = NULL, conserve = FALSE, conditioning = FALSE, .progressBar = TRUE) {

  b <- basisSet(modelList, direction)

  if(any(duplicated(names(b))) & conserve == FALSE & is.null(direction)) {

    dupOutput(b)

  }

  if(length(b) == 0) {

    warning("No independence claims present. Tests of directed separation not possible.", call. = FALSE)

    data.frame()

  } else {

    data <- modelList$data

    modelList <- removeData(modelList, formulas = 1)

    formulaList <- lapply(listFormula(modelList, formulas = 1), all.vars_trans)

    if(.progressBar == T & length(b) > 0)

      pb <- txtProgressBar(min = 0, max = length(b), style = 3)

    ret <- do.call(rbind, lapply(1:length(b), function(i) {

      bMod <- modelList[[which(sapply(formulaList, function(x) x[1] == b[[i]][2]))]]

      if(any(class(bMod) %in% c("lmerMod", "merModLmerTest", "glmerMod"))) {

        bNewMod <- suppressWarnings(
          update(bMod,
                 formula(paste(". ~ ", paste(rev(b[[i]][-2]), collapse = " + "), " + ", onlyBars(formula(bMod)))),
                 data = data)
        )

      } else

        bNewMod <- suppressWarnings(
          update(bMod,
                 formula(paste(". ~ ", paste(rev(b[[i]][-2]), collapse = " + "))),
                 data = data)
          )

      if(any(class(bNewMod) %in% c("lmerMod", "merModLmerTest"))) {

        kr <- KRp(bNewMod, b[[i]][1], data, intercepts = FALSE)

        ct <- summary(bNewMod)$coefficients

        ret <- data.frame(
          t(ct[which(b[[i]][1] == labels(terms(bNewMod))) + 1, 1:2]),
          kr[1, ],
          ct[which(b[[i]][1] == labels(terms(bNewMod))) + 1, 3],
          kr[2, ],
          row.names = NULL
        )

      }

      if(any(class(bNewMod) %in% c("lm", "glm", "negbin", "glmerMod"))) {

        ct <- as.data.frame(summary(bNewMod)$coefficients)

        ret <- ct[which(b[[i]][1] == labels(terms(bNewMod))) + 1, , drop = FALSE]

        if(all(class(bNewMod) %in% c("lm", "glm", "negbin"))) ret <- cbind(ret[, 1:2], DF = summary(bNewMod)$df[2], ret[, 3:4])

        if(all(class(bNewMod) %in% c("glmerMod", "pgls"))) ret <- cbind(ret[, 1:2], DF = NA, ret[, 3:4])

      }

      if(all(class(bNewMod) == "pgls")) {

        ct <- as.data.frame(summary(bNewMod)$coefficients)

        ret <- ct[which(b[[i]][1] == rownames(ct)), , drop = FALSE]

        ret <- cbind(ret[, 1:2], DF = bNewMod$n, ret[, 3:4])

      }

      if(any(class(bNewMod) %in% c("phylolm", "phyloglm"))) {

        ct <- as.data.frame(summary(bNewMod)$coefficients)

        ret <- ct[which(b[[i]][1] == rownames(ct)), , drop = FALSE]

        ret <- cbind(ret[, 1:2], DF = bNewMod$n, ret[, c(3, 6)])

      }

      if(all(class(bNewMod) %in% c("sarlm"))) {

        ct <- as.data.frame(summary(bNewMod)$Coef)

        ret <- ct[which(b[[i]][1] == rownames(ct)), , drop = FALSE]

        ret <- cbind(ret[, 1:2], DF = NA, ret[, 3:4])

      }

      if(any(class(bNewMod) %in% c("gls", "lme", "glmmPQL"))) {

        ct <- as.data.frame(summary(bNewMod)$tTable)

        ret <- ct[which(b[[i]][1] == labels(terms(bNewMod))) + 1, , drop = FALSE]

        if(ncol(ret) == 4 & all(class(bNewMod) %in% c("gls")))

          ret <- cbind(ret[, 1:2], DF = length(residuals(bNewMod)), ret[, 3:4])

      }

      names(ret) <- c("Estimate", "Std.Error", "DF", "Crit.Value", "P.Value")

      rhs <- paste0(b[[i]][-2], " ", collapse = " + ")

      if(conditioning == FALSE)

        rhs <- paste0(b[[i]][1], " + ...")

      ret <- data.frame(Independ.Claim = paste(b[[i]][2], " ~ ", rhs), ret)

      if(.progressBar == TRUE) setTxtProgressBar(pb, i)

      return(ret)

    } ) )

    if(.progressBar == TRUE) close(pb)

    rownames(ret) <- NULL

    if(conserve == TRUE) {

      ret = do.call(rbind, lapply(unique(names(b)), function(i) {

        r = ret[which(names(b) == i), ]

        r[which.min(r[, "P.Value"]), ]

      } ) )

    }

    ret <- cbind.data.frame(ret, sig = sapply(ret$P.Value, isSig))

    names(ret)[ncol(ret)] <- ""

    return(ret)

  }

}

dupOutput <- function(b, conserve = FALSE) {

  dup <- names(b)[which(duplicated(names(b)))]

  if(conserve == FALSE) {

    s <- paste("\nNon-linearities detected in the basis set where P-values are not symmetrical.",
               "\nThis can bias the outcome of the tests of directed separation.\n",

               "\nOffending independence claims:",

               lapply(dup, function(i) {

                 d <- b[names(b) %in% dup]

                 paste(
                   "\n", paste(d[[1]][2], "<-", d[[1]][1]), "*OR*",
                   paste(d[[1]][2], "->", d[[1]][1]), "\n"
                 )

               } ),

               "\nOption 1: Specify directionality using argument 'direction = c()' in 'summary'.\n",

               "\nOption 2: Remove path from the basis set by specifying as a correlated error using '%~~%' in 'psem'.\n",

               "\nOption 3 (recommended): Use argument 'conserve = TRUE' in 'summary' to compute both tests, and return the most conservative P-value.\n"

               )

    stop(s, call. = FALSE)

  }

}
