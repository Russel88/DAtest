% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DA.fri.R
\name{DA.fri}
\alias{DA.fri}
\title{Friedman Rank Sum test}
\usage{
DA.fri(
  data,
  predictor,
  paired = NULL,
  relative = TRUE,
  p.adj = "fdr",
  allResults = FALSE,
  ...
)
}
\arguments{
\item{data}{Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples}

\item{predictor}{The predictor of interest. Factor, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation}

\item{paired}{For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation}

\item{relative}{Logical. Should \code{data} be normalized to relative abundances. Default TRUE}

\item{p.adj}{Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details}

\item{allResults}{If TRUE will return raw results from the \code{friedman.test} function}

\item{...}{Additional arguments for the \code{friedman.test} function}
}
\value{
A data.frame with with results.
}
\description{
Apply friedman test to multiple features with one \code{predictor}
}
\examples{
# Creating random count_table, predictor, and paired variable
set.seed(4)
mat <- matrix(rnbinom(1500, size = 0.1, mu = 500), nrow = 100, ncol = 15)
rownames(mat) <- 1:100
pred <- c(rep("A", 5), rep("B", 5), rep("C", 5))
subject <- rep(1:5, 3)

# Running Friedman test on each feature
res <- DA.fri(data = mat, predictor = pred, paired = subject)
}
