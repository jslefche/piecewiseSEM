#' Retrieve tests of directed separation for a list of structural equations
#'
#' @param modelList a list of structural equations

dSep <- function(modelList, direction = NULL, conserve = FALSE, conditional = FALSE, .progressBar = TRUE) {

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

    formulaList <- lapply(listFormula(modelList, formulas = 1), all.vars.trans)

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

      } else if(all(class(bMod) %in% c("pgls"))) {

        bNewMod <- update(bMod,
                          formula(paste(". ~ ", paste(rev(b[[i]][-2]), collapse = " + "))),
                          data = data)

        } else bNewMod <- update(bMod,
                                 formula(paste(". ~ ", paste(rev(b[[i]][-2]), collapse = " + "))),
                                 data = data)

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

      if(conditional == FALSE)

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

               "\nOption 1: Specify directionality using argument 'direction = c()'.\n",

               "\nOption 2: Remove path from the basis set by specifying as a correlated error using '%~~%'.\n",

               "\nOption 3: Use argument 'conserve = TRUE' to compute both tests, and return the most conservative P-value.\n"

               )

    stop(s, call. = FALSE)

  }

}
