#' Retrieve tests of directed separation for a list of structural equations

#' @param modelList a list of structural equations

dSep <- function(modelList, conditional = FALSE, .progressBar = TRUE) {

  b <- basisSet(modelList)

  if(length(b) == 0)

    paste("No independence claims present. Tests of directed separation not possible.")

  else {

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

        kr <- KRmodcomp(bNewMod, bNewmod_drop)

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

      }

      if(any(class(bNewMod) %in% c("lme", "glmmPQL"))) {

        ct <- summary(bNewMod)$tTable

        ret <- ct[nrow(ct), ]

      }

      names(ret) <- c("Estimate", "Std.Error", "DF", "Crit.Value", "P.value")

      rhs <- paste(rev(b[[i]][-2]), collapse = " + ")

      if(conditional == FALSE)

        rhs <- paste(b[[i]][length(b[[i]])], " + ...")

      ret <- data.frame(Independence.Claim = paste(b[[i]][2], " ~ ", rhs), t(ret))

      if(.progressBar == TRUE) setTxtProgressBar(pb, i)

      return(ret)

    } ) )

    if(.progressBar == TRUE) close(pb)

    ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, is.numeric))], 4)

    ret <- cbind.data.frame(ret, sig = sapply(ret$P.value, isSig))

    names(ret)[ncol(ret)] <- ""

    return(ret)

  }

}
