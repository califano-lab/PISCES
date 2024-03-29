% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/narnea_funcs.R
\name{weighted_integration}
\alias{weighted_integration}
\title{Integrates 3D arrays of NaRnEA results generated by `narnea_combine` with sample-specific weights.
Weights are generated via a network voting procedure.}
\usage{
weighted_integration(combine.list, return.weights = TRUE)
}
\arguments{
\item{combine.list}{List generated by `narnea_combine`; two 3D arrays, one each for PES and NES}

\item{return.weights}{Flag to return network-sample weight matrix. TRUE by default.}
}
\value{
Integrated NaRnEA results, with PES and NES matrices.
}
\description{
Integrates 3D arrays of NaRnEA results generated by `narnea_combine` with sample-specific weights.
Weights are generated via a network voting procedure.
}
