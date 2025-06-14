% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{read_input_file}
\alias{read_input_file}
\title{Read input data from a file path}
\usage{
read_input_file(path)
}
\arguments{
\item{path}{Character string of the file path to read.}
}
\value{
The contents of the file, either a \code{data.frame}, \code{GRanges}, or other R object saved in the RDS.
}
\description{
Reads data from the specified file path. Supports RDS, CSV, TSV, and TXT formats.
Returns a data.frame or a GRanges object if RDS contains GRanges.
}
