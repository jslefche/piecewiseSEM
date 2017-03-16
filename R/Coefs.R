#' Get (standardized) coefficients from list of structural equations

#' @param modelList a list of structural equations

coefs <- function(modelList, data = NULL, intercepts = FALSE, standardize = TRUE) {

  if(class(modelList) == "psem") data <- modelList$data else

    data <- getData(modelList)

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)

  modelList <- modelList[!sapply(modelList, function(x) any(class(x) %in% c("matrix", "data.frame", "formula")))]

  tab <- getCoefs(modelList, data, intercepts)

  if(standardize == TRUE) tab <- data.frame(tab, Std.Estimate = stdCoefs(modelList, data, tab, intercepts))

  tab <- cbind(tab, isSig(tab$P.Value))

  names(tab)[length(tab)] = ""

  return(tab)

}

getCoefs <- function(modelList, data, intercepts = FALSE) {

  tab <- do.call(rbind, lapply(modelList, function(i) {

    if(all(class(i) %in% c("formula.cerror")))

      tab <- cerror(i, modelList, data) else {

      if(all(class(i) %in% c("lm", "glm", "lmerMod", "glmerMod", "merModLmerTest")))

        tab <- summary(i)$coefficients

      if(all(class(i) %in% c("lmerMod", "merModLmerTest"))) {

        krp <- KRp(i, all.vars.merMod(formula(i))[-1], intercepts = TRUE)

        tab <- as.data.frame(append(as.data.frame(tab), list(DF = krp[1,]), after = 2))

        tab[, "Pr(>|t|)"] <- krp[2, ]

        }

      if(all(class(i) %in% c("gls", "lme", "glmmPQL")))

        tab <- summary(i)$tTable

      tab <- data.frame(
        Response = all.vars.merMod(listFormula(list(i))[[1]])[1],
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

    } ) )

  tab[, which(sapply(tab, is.numeric))] <- round(tab[, which(sapply(tab, is.numeric))], 4)

  return(tab)

}

#' Get residual degrees of freedom for a linear regression
getDF <- function(model, tab) {

  cbind(tab[, 1:4], DF = NA, tab[, 5:6])

}

#' Calculate standardized regression coefficients
stdCoefs <- function(modelList, data, tab, intercepts) {

  newdata <- data

  do.call(c, lapply(1:length(modelList), function(i) {

    j <- modelList[[i]]

    f <- unlist(lapply(listFormula(list(j)), all.vars.merMod))

    if(all(class(j) %in% c("formula.cerror"))) {

      Bnew <- subset(tab, Response == paste0("~~", f[1]))$Estimate

    } else {

      notrans <- all.vars.notrans(formula(j))

      trans <- all.vars.trans(formula(j))

      if(any(notrans != trans)) {

        for(k in 1:length(notrans)) {

          newdata[, notrans[k]] <-

            sapply(newdata[, notrans[k]], function(x) eval(parse(text = gsub(notrans[k], x, trans[k]))))

        }

      }

      if(!identical(newdata, data)) tabNew <- getCoefs(j, newdata, intercepts) else tabNew <- tab

      B <- subset(tabNew, Response == f[1])$Estimate

      sd.x <- sapply(f[-1], function(x) sd(newdata[, x], na.rm = TRUE))

      sd.y <- sdFam(f[1], j, newdata)

      if(intercepts == FALSE) Bnew <- B * (sd.x / sd.y) else

        Bnew <- c(0, B[-1] * (sd.x / sd.y))

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

      y <- data[, x] else {

        if(any(class(model) %in% c("glmerMod")))

          linkfun <- model@resp$family$linkfun else

            linkfun <- model$family$linkfun

        y <- data[, x] * linkfun(mean(data[, x]))

      }

  }

  sd(y, na.rm = TRUE)

}
