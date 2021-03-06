% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meth_analysis.R
\name{dmTable}
\alias{dmTable}
\title{dmTable}
\usage{
dmTable(data, tested_samples, reference_samples)
}
\arguments{
\item{data}{A matrix of row methylation data with TCGA sample ids as column names and probes ids as row names.}

\item{tested_samples}{a vector of ids corresmasking to tumoral samples.}

\item{reference_samples}{a vector of ids corresmasking to control samples.}
}
\description{
Generate differential methylation data table from beta values. It basically extract all tumorous samples and controls. Then it computes the difference between each tumorous and the mean of control.
}

