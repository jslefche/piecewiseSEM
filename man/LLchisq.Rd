% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LLchisq.R
\name{LLchisq}
\alias{LLchisq}
\title{Generalized chi-squared for piecewise SEM}
\usage{
LLchisq(
  modelList,
  basis.set = NULL,
  direction = NULL,
  interactions = FALSE,
  conserve = FALSE
)
}
\arguments{
\item{modelList}{A list of structural equations created using \code{psem}.}

\item{basis.set}{An optional list of independence claims.}

\item{direction}{A \code{vector} of claims defining the specific
directionality of independence claims; for use in special cases (see
 \code{\link{dSep}}.}

\item{interactions}{whether interactions should be included in basis set. 
Default is FALSE}

\item{conserve}{Whether the most conservative log-likelihood should be returned;
for use in special cases (see Details). Default is FALSE.}
}
\value{
a data.frame corresponding to the Chi-squared statistic, d.f., and P-value
}
\description{
Derivation of log-likelihoods to be used in determining the
goodness-of-fit for piecewise structural equation models.
}
\details{
Here, a list of saturated models is first derived from the list of
structured equations using the basis set. Then, the differences in summed
log-likelihoods are computed and used to calculate the Chi-squared statistic.
}
\examples{
mod <- psem(
lm(rich ~ cover, data = keeley),
lm(cover ~ firesev, data = keeley),
lm(firesev ~ age, data = keeley),
data = keeley
)

LLchisq(mod)

}
\references{
Shipley, Bill, and Jacob C. Douma. "Generalized AIC and chi‐squared statistics 
for path models consistent with directed acyclic graphs." Ecology 101.3 (2020): e02960.
}
\seealso{
\code{\link{basisSet}}, \code{\link{dSep}}
}
\author{
Jon Lefcheck <LefcheckJ@si.edu>
}
