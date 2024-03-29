% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runtimeDA.R
\name{runtimeDA}
\alias{runtimeDA}
\title{Estimate runtime of \code{testDA} on large datasets}
\usage{
runtimeDA(
  data,
  predictor,
  paired = NULL,
  covars = NULL,
  subsamples = c(500, 1000, 1500, 2000),
  subsamples.slow = c(100, 150, 200, 250),
  tests = c("abc", "sam", "qua", "fri", "vli", "qpo", "pea", "wil", "ttt", "ttr",
    "ltt", "ltt2", "ere", "ere2", "msf", "zig", "lim", "lli", "lli2", "aov", "lao",
    "lao2", "kru", "lrm", "llm", "llm2", "spe", "aoa", "aoc", "tta", "ttc", "lma", "lmc",
    "lia", "lic"),
  tests.slow = c("mva", "neb", "bay", "per", "ds2", "ds2x", "zpo", "znb", "adx", "poi",
    "erq", "erq2"),
  cores = (detectCores() - 1),
  ...
)
}
\arguments{
\item{data}{Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples, and there should be rownames}

\item{predictor}{The predictor of interest. Either a Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation. If the \code{predictor} is numeric it will be treated as such in the analyses}

\item{paired}{For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation. Only for "poi", "per", "ttt", "ltt", "ltt2", "neb", "wil", "erq", "ds2", "lrm", "llm", "llm2", "lim", "lli", "lli2" and "zig"}

\item{covars}{Either a named list with covariates, OR if \code{data} is a \code{phyloseq} object a character vector with names of the variables in \code{sample_data(data)}}

\item{subsamples}{Vector with numbers of features to subsample to estimate runtime for fast methods}

\item{subsamples.slow}{Vector with numbers of features to subsample to estimate runtime for slow methods}

\item{tests}{Fast methods to include}

\item{tests.slow}{Slow methods to include}

\item{cores}{Integer. Number of cores to use for parallel computing. Default one less than available. Set to 1 for sequential computing.}

\item{...}{Additional arguments for the \code{testDA} function}
}
\value{
A data.frame with estimated runtimes for 1 run
}
\description{
Estimate the runtime of \code{testDA} from running on a subset of the features. Intended for datasets with at least 2000 features.
}
\details{
Outputs the estimated times for running each method 1 time. With cores=1 the runtime will be the sum of them all. With more cores the actual runtime will decrease asymptotically towards the slowest test

Runtime of all methods are expected to scale linearly with the number of features, except "anc" and "bay" which are modelled with a 2. order polynomial.
}
\examples{
# Creating large random count_table and predictor
set.seed(5)
mat <- matrix(rnbinom(150000, size = 0.5, mu = 500), nrow = 10000, ncol = 10)
rownames(mat) <- 1:10000
pred <- c(rep("A", 5), rep("B", 5))

# Use runtimeDA to predict total runtime for all features
# This example uses 1 core (cores = 1). 
# Remove the cores argument to get it as high (and thereby fast) as possible.
# Also, in this example only a subset of tests are run.
runtimeDA(mat, pred, cores = 1, tests = c("ttt","wil"), tests.slow = c("neb"))
}
