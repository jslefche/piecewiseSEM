#' #' Get (standardized) coefficients from list of structural equations
#'
#' #' @param modelList a list of structural equations
#'
#' Coefs <- function(modelList, data) {
#'
#'   tab <- lapply(modelList, summary)
#'
#' }
#'
#'
#' getCoefs <- function(modelList) {
#'
#'   do.call(rbind, lapply(modelList, function(i) {
#'
#'     if(all(class(i) %in% c("lm", "glm", "lmerMod", "glmerMod", "merModLmerTest"))) {
#'
#'       tab <- summary(i)$coefficients
#'
#'       if(all(class(i) %in% c("lmerMod", "merModLmerTest"))) {
#'
#'         KRp <- sapply(names(fixef(i))[-1], function(x) {
#'
#'           reducMod <- update(i, as.formula(paste("~ . -", x)))
#'
#'           pbkrtest::KRmodcomp(i, reducMod)$test$p.value[1]
#'
#'         } )
#'
#'         tab[, `Pr(>|t|)`] <- KRp
#'
#'       }
#'
#'     if(all(class(i) %in% c("gls", "nlme", "glmmPQL")))
#'
#'       tab <- summary(i)$tTable
#'
#'     return(tab)
#'
#'   } ) )
#'
#' }
#'
#'
#' stdCoefs <- function(modelList, data) {
#'
#'
#' }
