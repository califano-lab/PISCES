% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene-name_funcs.R
\name{GeneNameConvert}
\alias{GeneNameConvert}
\title{Converts gene names (rownames) from 'from' to 'to' convenctions.}
\usage{
GeneNameConvert(dat.mat, species, from, to)
}
\arguments{
\item{dat.mat}{Matrix with rownames to be changed.}

\item{species}{Species for this data. Either 'human' or 'mouse'.}

\item{from}{Starting name convention. One of 'ensg', 'gn', or 'entrez'.}

\item{to}{Ending name convention. One of 'ensg', 'gn', or 'entrez'.}
}
\value{
Matrix with converted row names.
}
\description{
Converts gene names (rownames) from 'from' to 'to' convenctions.
}
