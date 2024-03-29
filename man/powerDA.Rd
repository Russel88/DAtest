% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/powerDA.R
\name{powerDA}
\alias{powerDA}
\title{Estimating (empirical) statistical power}
\usage{
powerDA(
  data,
  predictor,
  paired = NULL,
  covars = NULL,
  test = NULL,
  effectSizes = c(2, 4, 8, 16, 32),
  alpha.p = 0.05,
  alpha.q = 0.1,
  p.adj = "fdr",
  R = 5,
  relative = TRUE,
  k = NULL,
  cores = (detectCores() - 1),
  args = list(),
  out.all = NULL,
  core.check = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{data}{Either a data.frame with counts/abundances, OR a \code{phyloseq} object. If a data.frame is provided rows should be taxa/genes/proteins and columns samples, and there should be rownames}

\item{predictor}{The predictor of interest. Either a Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation. If the \code{predictor} is numeric it will be treated as such in the analyses}

\item{paired}{For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation.}

\item{covars}{Either a named list with covariates, OR if \code{data} is a \code{phyloseq} object a character vector with names of the variables in \code{sample_data(data)}}

\item{test}{Character. Which test to include. See \code{testDA} for details on the implemented tests.}

\item{effectSizes}{Numeric. The effect sizes for the spike-ins. Default \code{c(2,4,8,16,32)}}

\item{alpha.p}{p-value threshold for false positive rates. Default 0.05}

\item{alpha.q}{q-value threshold for determining significance for \code{empirical power}. Default 0.1. This will change \code{fdr.output} for "sam".}

\item{p.adj}{Character. Method for p-value adjustment. See \code{p.adjust} for details. Default "fdr"}

\item{R}{Integer. Number of times to run the tests. Default 5}

\item{relative}{Logical. TRUE (default) for compositional data, FALSE for absolute abundances or pre-normalized data.}

\item{k}{Vector of length 3. Number of Features to spike in each tertile (lower, mid, upper). E.g. \code{k=c(5,10,15)}: 5 features spiked in low abundance tertile, 10 features spiked in mid abundance tertile and 15 features spiked in high abundance tertile. Default NULL, which will spike 2 percent of the total amount of features in each tertile (a total of 6 percent), but minimum c(5,5,5)}

\item{cores}{Integer. Number of cores to use for parallel computing. Default one less than available. Set to 1 for sequential computing.}

\item{args}{List. A list with arguments passed to method.}

\item{out.all}{If TRUE linear models will output results and p-values from \code{anova}/\code{drop1}, ds2/ds2x will run LRT and not Wald test, erq and erq2 will produce one p-value for the predictor, and limma will run F-tests. If FALSE will output results for 2. level of the \code{predictor}. If NULL (default) set as TRUE for multi-class predictors and FALSE otherwise}

\item{core.check}{If TRUE will make an interactive check that the amount of cores specified are desired. Only if \code{cores>20}. This is to ensure that the function doesn't automatically overloads a server with workers.}

\item{verbose}{If TRUE will print informative messages}
}
\value{
An object of class \code{DAPower}, which contains a list with 1: A data.frame with results, 2: alpha.p value, 3: alpha.q value
}
\description{
Estimating (empirical) statistical power for a specific differential abundance and expression method on a specific dataset
}
\details{
Currently implemented methods: see \code{testDA}
}
\examples{
# Creating random count_table and predictor
set.seed(5)
mat <- matrix(rnbinom(1000, size = 0.5, mu = 500), nrow = 50, ncol = 20)
rownames(mat) <- 1:50
pred <- c(rep("Control", 10), rep("Treatment", 10))

# Running powerDA on Wilcoxon test to test it with different effect sizes
# This example uses 1 core (cores = 1). 
# Remove the cores argument to get it as high (and thereby fast) as possible.
res <- powerDA(data = mat, predictor = pred, test = "wil", cores = 1)
summary(res)

\donttest{
# Include a paired variable for dependent/blocked samples
subject <- rep(1:10, 2)
res <- powerDA(data = mat, predictor = pred, paired = subject, test = "ttt", cores = 1)

# Include covariates
covar1 <- rnorm(20)
covar2 <- rep(c("A","B"), 10)
res <- powerDA(data = mat, predictor = pred, 
               covars = list(FirstCovar = covar1, CallItWhatYouWant = covar2), 
               test = "lrm", cores = 1)

# Data is absolute abundance
res <- powerDA(data = mat, predictor = pred, relative = FALSE, test = "ttt", cores = 1)
}
}
