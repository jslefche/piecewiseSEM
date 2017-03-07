#' Retrieve tests of directed separation for a list of structural equations

#' @param modelList a list of structural equations

dSep <- function(modelList, direction = NULL, conserve = FALSE, conditional = FALSE, .progressBar = TRUE) {

  b <- basisSet(modelList, direction)

  if(length(b) == 0) {

    warning("No independence claims present. Tests of directed separation not possible.", call. = FALSE)

    data.frame()

  }

  else if(any(duplicated(names(b))) & conserve == FALSE & is.null(direction)) {

    dupOutput(b)

  } else {

    formulaList <- lapply(listFormula(modelList, remove = TRUE), all.vars.merMod)

    if(.progressBar == T & length(b) > 0)

      pb <- txtProgressBar(min = 0, max = length(b), style = 3)

    ret <- do.call(rbind, lapply(1:length(b), function(i) {

      bMod <- modelList[[which(sapply(formulaList, function(x) x[1] == b[[i]][2]))]]

      if(any(class(bMod) %in% c("lmerMod", "merModLmerTest", "glmerMod"))) {

        bNewMod <- suppressWarnings(
          update(bMod, formula(paste(". ~ ", paste(rev(b[[i]][-2]), collapse = " + "), " + ", onlyBars(formula(bMod)))))
        )

      } else {

        bNewMod <- update(bMod, formula(paste(". ~ ", paste(rev(b[[i]][-2]), collapse = " + "))))

      }

      if(any(class(bNewMod) %in% c("lmerMod", "merModLmerTest"))) {

        bNewmod_drop <- update(bNewMod, formula(paste(" ~ . - ", b[[i]][1])))

        kr <- pbkrtest::KRmodcomp(bNewMod, bNewmod_drop)

        ct <- summary(bNewMod)$coefficients

        ret <- data.frame(
          t(ct[nrow(ct), 1:2]),
          kr$test$ddf[1],
          ct[nrow(ct), 3],
          kr$test$p.value[1],
          row.names = NULL
        )

      }

      if(any(class(bNewMod) %in% c("lm", "glm", "glmerMod"))) {

        ct <- summary(bNewMod)$coefficients

        ret <- ct[nrow(ct), ]

        ret <- c(ret[1:2], DF = NA, ret[3:4])

        # Add in df

      }

      if(any(class(bNewMod) %in% c("lme", "glmmPQL"))) {

        ct <- summary(bNewMod)$tTable

        ret <- ct[nrow(ct), ]

      }

      names(ret) <- c("Estimate", "Std.Error", "DF", "Crit.Value", "P.value")

      rhs <- paste0(b[[i]][-2], " ", collapse = "+")

      if(conditional == FALSE)

        rhs <- paste0(b[[i]][1], " +...")

      ret <- data.frame(Independ.Claim = paste(b[[i]][2], " ~ ", rhs), t(ret))

      if(.progressBar == TRUE) setTxtProgressBar(pb, i)

      return(ret)

    } ) )

    if(.progressBar == TRUE) close(pb)

    if(conserve == TRUE) {

      ret = do.call(rbind, lapply(unique(names(b)), function(i) {

        r = ret[which(names(b) == i), ]

        r[which.min(r[, "P.value"]), ]

      } ) )

    }

    ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, is.numeric))], 4)

    ret <- cbind.data.frame(ret, sig = sapply(ret$P.value, isSig))

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
