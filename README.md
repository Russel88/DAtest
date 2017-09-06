DAtest
======

This is a package for comparing different differential abundance methods
used in microbial marker-gene (e.g. 16S rRNA), RNA-seq and protein
abundance analysis.

There are many methods for testing differential abundance and no golden
standard, but results can vary a lot between the different statistical
methods. The false positive rate and the power of the different methods
highly depends on the dataset. This package aims at aiding the analyst
in choosing a method based on empirical testing.

#### The method goes as follows:

-   Shuffle outcome variable (E.g. case vs. control)
-   Spike in data for some randomly chosen features (OTU/gene/protein),
    such that they are associated with the shuffled outcome
-   Apply methods, and check:
    -   whether they can find the spike-ins
    -   whether the false positive rate is controlled

#### The workflow (details can be found below):

-   Compare methods with testDA function
    -   Check the results with `plot` or `summary`
    -   Choose method that has high AUC, and FPR not higher than ~0.05
-   Run data with the chosen test with DA."test" function, where "test"
    is the name of the test (see details with ?testDA)
    -   Check out your final results.

### Citation

Some scripts, including the spike-in for estimating AUC, is borrowed
from: [Thorsen J, Brejnrod A et al. Large-scale benchmarking reveals
false discoveries and count transformation sensitivity in 16S rRNA gene
amplicon data analysis methods used in microbiome studies. *Microbiome*
(2016)](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8)

Overview of this tutorial
-------------------------

-   [Installation of packages](#installation-of-packages)
-   [How to compare methods](#how-to-compare-methods)
-   [How to run real (unshuffled)
    data](#how-to-run-real-(unshuffled)-data)
-   [Implemented methods](#implemented-methods)
-   [Extra features](#extra-features)
-   [Errors and issues](#errors-and-issues)

Installation of packages
========================

    library(devtools)
    install_github("Russel88/DAtest")

#### The following are needed for *full* functionality

But the package will work without them

    source("https://bioconductor.org/biocLite.R")
    biocLite("DESeq2")
    biocLite("edgeR")
    biocLite("metagenomeSeq")
    biocLite("baySeq")
    biocLite("ALDEx2")
    biocLite("limma")

-   RAIDA has to be installed from an
    [external source.](https://cals.arizona.edu/~anling/software/)
    -   It depends on: MASS, protoclust, qvalue (biocLite("qvalue")) and
        limma (biocLite("limma"))

How to compare methods
======================

First, all methods with a "False Positive Rate" (FPR) above ~0.05 has an
inflated false positive rate, and the p-values can therefore not be
trusted.

Second, among methods with a low FPR, those with a high "Area Under the
(Receiver Operator) Curve" (AUC) and "Spike Detection Rate"
(Spike.detect.rate) are likely to have more power to detect signal in
the respective dataset.

Therefore, we want a method with a FPR ~0.05 or lower and a high AUC.

### **Run the test:**

(if you have a phyloseq object, scroll down)

    mytest <- testDA(data, outcome)

data is either a matrix or data.frame with taxa/genes/proteins as rows
and samples as columns.

outcome is the outcome of interest, e.g. a factor denoting whether
samples are cases or controls (in the same order as columns in data).

outcome can be a factor with more than two levels, in which case only
the second level is spiked, and if methods output several p-values, only
the p-value associated with the second level is used.

outcome can also be numeric.

**The function automatically uses multiple CPUs for fast execution**

The methods run in parallel, and by default the number of cores used is
one less than available. This can be changed with the `cores` argument.
`cores = 1` will turn off parallel computing. If the function is
terminated before ending you might get the following warning:

    closing unused connection X (...)

This can safely be ignored.

If you have terminated the function before it ended and your computer
runs slow, you might wan't to call `closeAllConnections()`. **WARNING:**
This might affect other R process you have running! Alternatively,
restart R.

### *If you have a paired/blocked experimental design:*

E.g. repeated samples from same patients, or control/treatments inside
blocks.

The paired argument should be a factor with the ID of the
patient/subject/block (in the same order as columns in data):

    mytest <- testDA(data, outcome, paired = SubjectID)

When a *paired* argument is provided, the outcome is shuffled within the
levels of the *paired* factor.

Paired analysis can be very slow. If you simply can't wait to see the
results remove "neb" from the tests argument.

### *If you have non-relative abundances, e.g. for normalized protein abundance:*

    mytest <- testDA(data, outcome, relative = FALSE)

#### *If you have negative values in your abundance table:*

    mytest <- testDA(data, outcome, relative = FALSE, spikeMethod = "add")

It might be an advantage to change the test statistic of the permutation
test for non-relative data, and this is highly recommended if the data
contains negative values as the default log ratio then becomes
nonsensical.

Example (alternatively use medians, if data is highly non-normal):

    testFun <- function(case,control){
      mean(control) - mean(case)
    }

    testFun.pair <- function(case,control){
      mean(control - case)
    }

    mytest <- testDA(data,outcome,relative=FALSE, 
                     args = list(per = list(testStat=testFun, testStat.pair=testFun.pair)))

Note that the output from ttt, ltt, ltt2 and wil also contains the log
ratio fold change. This does not affect the results of the testDA
function, but if any of these are used in the final analysis it should
be changed if the data contains negative values.

### If you have a phyloseq object:

data can also be a phyloseq object. In this case, the outcome and paired
arguments are the names of the variables in sample\_data(data):

    mytest <- testDA(data, outcome = "Time", paired = "Patient")

**Plot the output:**
--------------------

    plot(mytest, sort = "AUC")

**Print the output:**
---------------------

Medians for each method:

    summary(mytest, sort = "AUC")

How to run real (unshuffled) data
=================================

All tests can easily be run with the original data. E.g. edgeR exact
test:

    res.ere <- DA.ere(data, outcome)

All methods can be accessed in the same way; DA."test" where "test" is
the abbreviation given in the details of the testDA function.

Alternatively, run all (or several) methods and check which features are
found by several methods

    res.all <- allDA(data, outcome)

    head(res.all$table)

Implemented methods
===================

-   per - [Permutation test with user defined test
    statistic](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8)
    (See below for description of the paired permutation test)
-   bay -
    [baySeq](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-422)
-   adx - [ALDEx t-test and
    wilcoxon](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0067019)
-   wil - Wilcoxon Rank Sum on relative abundances (Paired is a Wilcoxon
    Signed Rank test)
-   ttt - Welch t.test on relative abundances
-   ltt - Welch t.test, but reads are first transformed with
    log(abundance + delta) then turned into relative abundances
-   ltt2 - Welch t.test, but with relative abundances transformed with
    log(relative abundance + delta)
-   neb - Negative binomial GLM with log of library size as offset (The
    paired version is a mixed-effect model)
-   erq - [EdgeR - Quasi
    likelihood](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/)
    (The paired version is a model with the paired variable
    as covariate)
-   ere - [EdgeR - Exact
    test](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/)
-   msf - [MetagenomeSeq feature
    model](https://bioconductor.org/packages/release/bioc/html/metagenomeSeq.html)
-   zig - [MetagenomeSeq zero-inflated
    gaussian](https://www.nature.com/nmeth/journal/v10/n12/full/nmeth.2658.html)
-   ds2 -
    [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)
    (The paired version is a model with the paired variable
    as covariate)
-   lim -
    [LIMMA](https://link.springer.com/chapter/10.1007%2F0-387-29362-0_23?LI=true)
    (The paired version is a model with the paired variable
    as covariate)
-   lli -
    [LIMMA](https://link.springer.com/chapter/10.1007%2F0-387-29362-0_23?LI=true),
    but reads are first transformed with log(abundance + delta) then
    turned into relative abundances
-   lli2 -
    [LIMMA](https://link.springer.com/chapter/10.1007%2F0-387-29362-0_23?LI=true),
    but with relative abundances transformed with log(relative
    abundance + delta)
-   kru - Kruskal-Wallis test on relative abundances
-   aov - ANOVA on relative abundances
-   lao - ANOVA, but reads are first transformed with log(abundance +
    delta) then turned into relative abundances
-   lao2 - ANOVA, but with relative abundances transformed with
    log(relative abundance + delta)
-   lrm - Linear regression on relative abundances (The paired version
    is a mixed-effect model)
-   llm - Linear regression, but reads are first transformed with
    log(abundance + delta) then turned into relative abundances
-   llm2 - Linear regression, but with relative abundances transformed
    with log(relative abundance + delta)
-   rai -
    [RAIDA](https://academic.oup.com/bioinformatics/article/31/14/2269/256302/A-robust-approach-for-identifying-differentially?searchresult=1)
-   spe - Spearman Rank Correlation

### Paired permutation test

A paired permutation test is implemented specifically for this package.
The test is similar to the original, but with a different test statistic
and permutation scheme. The permutations are constrained in the paired
version such that the outcome is only permuted within each level of the
paired argument (e.g. subjects). The test statistic first finds the
log-ratio between the two outcome levels (e.g. case and control) for
each level of the paired argument and the final statistic is the mean of
these log-ratios.

Extra features
==============

#### Plot the p-value distributions. Raw p-values should in theory have a uniform (flat) distribution between 0 and 1.

    plot(mytest, p = TRUE)

#### Results from all the runs:

    print(mytest)

#### See the output from the individual methods. E.g. "ere" first run:

    View(mytest$results[[1]]["ere"])

### Passing arguments to the different tests

Additional arguments can be passed to the internal functions with the
"args" argument. It should be structured as a list with elements named
by the tests:

E.g. passing to the DA.per function that it should only run 1000
iterations:

    mytest <- testDA(...,args = list(per=list(noOfIterations=1000)))

Include that the log t.test should use a pseudocount of 0.1:

    mytest <- testDA(...,args = list(per=list(noOfIterations=1000), ltt=list(delta=0.1)))

Additional arguments are simply seperated by commas.

Below is an overview of which functions get the arguments that are
passed to a specific test:

-   per - Passed to DA.per
-   bay - Passed to getPriors and getLikelihoods
-   adx - Passed to aldex
-   wil - Passed to wilcox.test and DA.wil
-   ttt - Passed to t.test and DA.ttt
-   ltt - Passed to t.test and DA.ltt
-   ltt2 - Passed to t.test and DA.ltt2
-   neb - Passed to glm.nb and glmer.nb
-   erq - Passed to exactTest
-   ere - Passed to glmQLFit
-   msf - Passed to fitFeatureModel
-   zig - Passed to fitZig
-   ds2 - Passed to DESeq
-   lim - Passed to eBayes
-   lli - Passed to eBayes
-   lli2 - Passed to eBayes
-   kru - Passed to kruskal.test
-   aov - Passed to aov
-   lao - Passed to aov
-   lao2 - Passed to aov
-   lrm - Passed to lm and lme
-   llm - Passed to lm and lme
-   llm2 - Passed to lm and lme
-   rai - Passed to raida
-   spe - Passed to cor.test

Errors and issues
=================

If a method fails the following is usually printed: *Error in { : task 1
failed - "task X failed"*. To find the method corresponding to X, run
the function again with *verbose = TRUE*. This will print the order of
the tests, and test number X can then be excluded from the *tests*
argument.
