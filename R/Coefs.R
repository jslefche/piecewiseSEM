#' Get (standardized) coefficients from list of structural equations

#' @param modelList a list of structural equations

coefs <- function(modelList, data) {

  getCoefs(modelList)

}

getCoefs <- function(modelList) {

  tab <- do.call(rbind, lapply(modelList, function(i) {

    if(all(class(i) %in% c("formula.cerror")))

      tab <- cerror(i, modelList) else {

      if(all(class(i) %in% c("lm", "glm", "lmerMod", "glmerMod", "merModLmerTest")))

        tab <- summary(i)$coefficients

      if(all(class(i) %in% c("lmerMod", "merModLmerTest"))) {

        KRp <- sapply(names(fixef(i))[-1], function(x) {

          reducMod <- update(i, as.formula(paste("~ . -", x)))

          pbkrtest::KRmodcomp(i, reducMod)$test$p.value[1]

        } )

        tab[, `Pr(>|t|)`] <- KRp

        }

      if(all(class(i) %in% c("gls", "nlme", "glmmPQL")))

        tab <- summary(i)$tTable

      tab <- data.frame(
        Response = unlist(lapply(listFormula(list(i)), all.vars.merMod))[1],
        Predictor = rownames(tab),
        tab,
        row.names = NULL
      )

      if(ncol(tab) == 6) {

        tab <- getDF(i, tab)

      }

      names(tab) <- c("Response", "Predictor", "Estimate", "Std.Error", "DF", "Crit.Value", "P.Value")

      }

    return(tab)

    }

  ) )

  tab[, which(sapply(tab, is.numeric))] <- round(tab[, which(sapply(tab, is.numeric))], 4)

  return(tab)

}

#' Get residual degrees of freedom for a linear regression
getDF <- function(model, tab) {

  cbind(tab[, 1:4], DF = NA, tab[, 5:6])

}


stdCoefs <- function(modelList, data) {


}
