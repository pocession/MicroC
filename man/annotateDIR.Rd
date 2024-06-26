% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotateDIR.R
\name{annotateDIR}
\alias{annotateDIR}
\title{annotateDIR}
\usage{
annotateDIR(input, output)
}
\arguments{
\item{input}{a character string specifying the name and path of the input file (.csv)
This file should be generated by getDIRWithNoReplicate()}

\item{output}{A character string specifying the name and path of output file
the default path is as same as the treatment file
the default output file name is the treatment file name prefixed with "annotated_"}
}
\value{
a dataframe of the output data
}
\description{
This function reads differentially-interacted regions (DIRs) and annotate with hg38 information
}
\examples{
\dontrun{
df <- annotateDIR(
input = here::here("./Results/DIR/bg.csv"),
output = NULL
)
}

}
\author{
Tsunghan Hsieh
}
