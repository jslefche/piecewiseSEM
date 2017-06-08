\name{rsquared}
\alias{rsquared}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rsquared(modelList, method = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelList}{
%%     ~~Describe \code{modelList} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
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
function (modelList, method = NULL) 
{
    if (!all(class(modelList) \%in\% c("psem", "list"))) 
        modelList = list(modelList)
    modelList <- removeData(modelList, formulas = 1)
    evaluateClasses(modelList)
    ret <- lapply(modelList, function(i) {
        if (any(class(i) \%in\% c("lm", "pgls"))) 
            r <- rsquared.lm(i)
        if (all(class(i) \%in\% c("gls"))) 
            r <- rsquared.gls(i)
        if (all(class(i) \%in\% c("sarlm"))) 
            r <- list(family = "gaussian", identity = "link", 
                R.squared = NA)
        if (any(class(i) \%in\% c("glm"))) 
            r <- rsquared.glm(i, method)
        if (all(class(i) \%in\% c("lme"))) 
            r <- rsquared.lme(i)
        if (all(class(i) \%in\% c("lmerMod", "merModLmerTest"))) 
            r <- rsquared.merMod(i)
        if (any(class(i) \%in\% c("glmerMod"))) 
            r <- rsquared.glmerMod(i, method)
        if (any(class(i) \%in\% c("glmmPQL"))) 
            r <- rsquared.glmmPQL(i, method)
        ret <- do.call(data.frame, r)
        ret <- data.frame(Response = all.vars.merMod(formula(i))[1], 
            ret)
        return(ret)
    })
    if (length(unique(sapply(ret, ncol))) != 1) {
        nc <- max(sapply(ret, ncol))
        ret <- lapply(ret, function(i) if (ncol(i) == nc) 
            i
        else {
            data.frame(i[, 1:3], Marginal = i[, 4], Conditional = NA)
        })
    }
    ret <- do.call(rbind, ret)
    rownames(ret) <- NULL
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
