\name{summary.psem}
\alias{summary.psem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.psem(modelList, groups = NULL, direction = NULL, conserve = FALSE, conditional = FALSE, add.claims = NULL, intercepts = FALSE, standardize = TRUE, .progressBar = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelList}{
%%     ~~Describe \code{modelList} here~~
}
  \item{groups}{
%%     ~~Describe \code{groups} here~~
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
  \item{add.claims}{
%%     ~~Describe \code{add.claims} here~~
}
  \item{intercepts}{
%%     ~~Describe \code{intercepts} here~~
}
  \item{standardize}{
%%     ~~Describe \code{standardize} here~~
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
function (modelList, groups = NULL, direction = NULL, conserve = FALSE, 
    conditional = FALSE, add.claims = NULL, intercepts = FALSE, 
    standardize = TRUE, .progressBar = TRUE) 
{
    if (!is.null(groups)) {
        data <- modelList$data
        lvls <- unique(data[, colnames(data) \%in\% groups])
        modelList2 <- lapply(lvls, function(i) update(modelList, 
            data = data[data[, colnames(data) \%in\% groups] == 
                i, ]))
        names(modelList2) <- lvls
        lapply(modelList2, function(i) summary.psem2(i, direction, 
            conserve, conditional, add.claims, intercepts, standardize, 
            .progressBar))
    }
    else summary.psem2(modelList, direction, conserve, conditional, 
        add.claims, intercepts, standardize, .progressBar)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
