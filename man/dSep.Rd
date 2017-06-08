\name{dSep}
\alias{dSep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dSep(modelList, direction = NULL, conserve = FALSE, conditional = FALSE, .progressBar = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelList}{
%%     ~~Describe \code{modelList} here~~
}
  \item{direction}{
%%     ~~Describe \code{direction} here~~
}
  \item{conserve}{
%%     ~~Describe \code{conserve} here~~
}
  \item{conditional}{
%%     ~~Describe \code{conditional} here~~
}
  \item{.progressBar}{
%%     ~~Describe \code{.progressBar} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (modelList, direction = NULL, conserve = FALSE, conditional = FALSE, 
    .progressBar = TRUE) 
{
    b <- basisSet(modelList, direction)
    if (any(duplicated(names(b))) & conserve == FALSE & is.null(direction)) {
        dupOutput(b)
    }
    if (length(b) == 0) {
        warning("No independence claims present. Tests of directed separation not possible.", 
            call. = FALSE)
        data.frame()
    }
    else {
        data <- modelList$data
        formulaList <- lapply(listFormula(modelList, remove = TRUE), 
            all.vars.trans)
        if (.progressBar == T & length(b) > 0) 
            pb <- txtProgressBar(min = 0, max = length(b), style = 3)
        ret <- do.call(rbind, lapply(1:length(b), function(i) {
            bMod <- modelList[[which(sapply(formulaList, function(x) x[1] == 
                b[[i]][2]))]]
            if (any(class(bMod) \%in\% c("lmerMod", "merModLmerTest", 
                "glmerMod"))) {
                bNewMod <- suppressWarnings(update(bMod, formula(paste(". ~ ", 
                  paste(rev(b[[i]][-2]), collapse = " + "), " + ", 
                  onlyBars(formula(bMod)))), data = data))
            }
            else if (all(class(bMod) \%in\% c("pgls"))) {
                bNewMod <- update(bMod, formula(paste(". ~ ", 
                  paste(rev(b[[i]][-2]), collapse = " + "))), 
                  data = data)
            }
            else bNewMod <- update(bMod, formula(paste(". ~ ", 
                paste(rev(b[[i]][-2]), collapse = " + "))), data = data)
            if (any(class(bNewMod) \%in\% c("lmerMod", "merModLmerTest"))) {
                kr <- KRp(bNewMod, b[[i]][1], data, intercepts = FALSE)
                ct <- summary(bNewMod)$coefficients
                ret <- data.frame(t(ct[nrow(ct), 1:2]), kr[1, 
                  ], ct[nrow(ct), 3], kr[2, ], row.names = NULL)
            }
            if (any(class(bNewMod) \%in\% c("lm", "glm", "glmerMod", 
                "pgls"))) {
                ct <- as.data.frame(summary(bNewMod)$coefficients)
                ret <- ct[nrow(ct), , drop = FALSE]
                if (all(class(bNewMod) \%in\% c("lm", "glm"))) 
                  ret <- cbind(ret[, 1:2], DF = summary(bNewMod)$df[2], 
                    ret[, 3:4])
                if (all(class(bNewMod) \%in\% c("glmerMod", "pgls"))) 
                  ret <- cbind(ret[, 1:2], DF = NA, ret[, 3:4])
            }
            if (all(class(bNewMod) \%in\% c("sarlm"))) {
                ct <- as.data.frame(summary(bNewMod)$Coef)
                ret <- ct[nrow(ct), , drop = FALSE]
                ret <- cbind(ret[, 1:2], DF = NA, ret[, 3:4])
            }
            if (any(class(bNewMod) \%in\% c("gls", "lme", "glmmPQL"))) {
                ct <- as.data.frame(summary(bNewMod)$tTable)
                ret <- ct[nrow(ct), , drop = FALSE]
                if (ncol(ret) == 4 & class(bNewMod) \%in\% c("gls")) 
                  ret <- cbind(ret[, 1:2], DF = length(residuals(bNewMod)), 
                    ret[, 3:4])
                if (ncol(ret) == 4 & class(bNewMod) \%in\% c("lme", 
                  "glmmPQL")) 
                  ret <- cbind(ret[, 1:2], DF = summary(bNewMod)$fixDF$X[-1], 
                    ret[, 3:4])
            }
            names(ret) <- c("Estimate", "Std.Error", "DF", "Crit.Value", 
                "P.Value")
            rhs <- paste0(b[[i]][-2], " ", collapse = " + ")
            if (conditional == FALSE) 
                rhs <- paste0(b[[i]][1], " + ...")
            ret <- data.frame(Independ.Claim = paste(b[[i]][2], 
                " ~ ", rhs), ret)
            if (.progressBar == TRUE) 
                setTxtProgressBar(pb, i)
            return(ret)
        }))
        if (.progressBar == TRUE) 
            close(pb)
        rownames(ret) <- NULL
        if (conserve == TRUE) {
            ret = do.call(rbind, lapply(unique(names(b)), function(i) {
                r = ret[which(names(b) == i), ]
                r[which.min(r[, "P.Value"]), ]
            }))
        }
        ret <- cbind.data.frame(ret, sig = sapply(ret$P.Value, 
            isSig))
        names(ret)[ncol(ret)] <- ""
        return(ret)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
