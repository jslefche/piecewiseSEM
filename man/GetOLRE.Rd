% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{GetOLRE}
\alias{GetOLRE}
\title{Obtain (observation-level) random effects from a generalized linear mixed model}
\usage{
GetOLRE(sigma, model, X, data, RE = c("all", "RE", "OLRE"))
}
\description{
RE = "all" all random effects are reported
RE = "RE" just group effects are reported
RE = "OLRE" just observation-level effects are reported
}
\keyword{internal}
