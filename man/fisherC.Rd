% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fisherC.R
\name{fisherC}
\alias{fisherC}
\title{Summarize tests of directed separation using Fisher's C statistic}
\usage{
fisherC(
  dTable,
  add.claims = NULL,
  basis.set = NULL,
  direction = NULL,
  interactions = FALSE,
  conserve = FALSE,
  conditional = FALSE,
  .progressBar = FALSE
)
}
\arguments{
\item{dTable}{a \code{data.frame} containing tests of directed separation from \code{dSep}}

\item{add.claims}{an optional vector of additional independence claims (i.e., P-values) 
to be added to the basis set}

\item{basis.set}{An optional list of independence claims.}

\item{direction}{a vector of claims defining the specific directionality of any independence 
claim(s)}

\item{interactions}{whether interactions should be included in independence claims. 
Default is FALSE}

\item{conserve}{whether the most conservative P-value should be returned. 
Default is FALSE}

\item{conditional}{whether the conditioning variables should be shown in the table. 
Default is FALSE}

\item{.progressBar}{an optional progress bar. Default is FALSE}
}
\value{
a data.frame corresponding to the C statistic, d.f., and P-value
}
\description{
Summarize tests of directed separation using Fisher's C statistic
}
