#' Get (standardized) coefficients from list of structural equations

#' @param modelList a list of structural equations

coefs <- function(modelList, data, intercepts = FALSE, standardize = TRUE) {

  tab <- getCoefs(modelList, data, intercepts)

  if(standardize == TRUE) tab <- data.frame(tab, Std.Estimate = stdCoefs(modelList, tab, data, intercepts))

  tab <- cbind(tab, isSig(tab$P.Value))

  names(tab)[length(tab)] = ""

  return(tab)

}

getCoefs <- function(modelList, data, intercepts = FALSE) {

  modelList <- modelList[!sapply(modelList, function(x) any(class(x) == "formula"))]

  tab <- do.call(rbind, lapply(modelList, function(i) {

    if(all(class(i) %in% c("formula.cerror")))

      tab <- cerror(i, modelList, data) else {

      if(all(class(i) %in% c("lm", "glm", "lmerMod", "glmerMod", "merModLmerTest")))

        tab <- summary(i)$coefficients

      if(all(class(i) %in% c("lmerMod", "merModLmerTest"))) {

        KRp <- sapply(names(fixef(i))[-1], function(x) {

          reducMod <- update(i, as.formula(paste("~ . -", x)))

          pbkrtest::KRmodcomp(i, reducMod)$test$p.value[1]

        } )

        tab[, `Pr(>|t|)`] <- KRp

        }

      if(all(class(i) %in% c("gls", "lme", "glmmPQL")))

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

    if(intercepts == FALSE) tab <- subset(tab, Predictor != "(Intercept)")

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

#' Calculate standardized regression coefficients
stdCoefs <- function(modelList, tab = NULL, data, intercepts) {

  modelList <- modelList[!sapply(modelList, function(x) any(class(x) == "formula"))]

  if(is.null(tab)) tab <- getCoefs(modelList, data, intercepts)

  do.call(c, lapply(1:length(modelList), function(i) {

    if(is.list(data)) newdata <- data[[i]] else newdata <- data

    i <- modelList[[i]]

    notrans <- all.vars.notrans(formula(i))

    trans <- all.vars.trans(formula(i))

    if(any(notrans != trans)) {

      for(j in 1:length(notrans)) {

        newdata[, notrans[j]] <-

          sapply(newdata[, notrans[j]], function(x) eval(parse(text = gsub(notrans[j], x, trans[j]))))

      }

    }

    f <- unlist(lapply(listFormula(list(i)), all.vars.merMod))

    if(all(class(i) %in% c("formula.cerror"))) {

      Bnew <- subset(tab, Response == f[1])$Estimate

    } else {

      B <- subset(tab, Response == f[1])$Estimate

      sd.x <- sapply(f[-1], function(x) sd(newdata[, x], na.rm = TRUE))

      sd.y <- sdFam(f[1], i, newdata)

      Bnew <- B * (sd.x / sd.y)

    }

    Bnew <- round(Bnew, 4)

    unname(Bnew)

  } ) )

}

#' Properly scale standard deviations depending on the error distribution
sdFam <- function(x, model, data) {

  .family = try(family(model), silent = TRUE)

  if(class(.family) == "try-error") y <- data[, x] else {

    if(.family$family != "gaussian")

      warning(
        paste0("Standardized coefficients for non-normal distributions must be interpreted differently. See help('stdCoefs')."),
        call. = FALSE
        )

    .link <- .family$link

    if(.link == "identity")

      y <- data[, x] else

        y <- data[, x] * model$family$linkfun(mean(data[, x]))

  }

  sd(y, na.rm = TRUE)

}
