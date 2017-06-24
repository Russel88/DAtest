DAtest
======

This is a package for comparing different differential abundance methods
used in both microbiome and RNA-seq analysis.

Most scripts, including the spike-in for estimating AUC, is borrowed
from: [Thorsen J, Brejnrod A et al. Large-scale benchmarking reveals
false discoveries and count transformation sensitivity in 16S rRNA gene
amplicon data analysis methods used in microbiome studies. *Microbiome*
(2016)](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8)

Installation of packages
------------------------

    library(devtools)
    install_github("Russel88/DAtest")

#### The following are needed for full functionality

    source("https://bioconductor.org/biocLite.R")
    biocLite("DESeq2")
    biocLite("edgeR")
    biocLite("metagenomeSeq")
    biocLite("baySeq")
    biocLite("ALDEx2")

ANCOM has to be installed from an [external
source,](https://www.niehs.nih.gov/research/resources/software/biostatistics/ancom/index.cfm)
and is not included by default in the test.

How to compare methods:
-----------------------

A good method has a "False Positive Rate" (FPR) at 0.05 or below and an
"Area Under the Curve" (AUC) as high as possible.

**Run the test:**

    mytest <- testDA(otu_table,predictor)

otu\_table is an OTU table with taxa as rows.

predictor is the outcome of interest, e.g. a factor denoting whether
samples are cases or controls (in the same order as columns in
otu\_table).

**Print the output:**

    summary(mytest)

**Plot the output:**

    plot(mytest)

How to run real data:
---------------------

All tests can easily be run with the original data. E.g. edgeR exact
test:

    res.ere <- DA.ere(otu_table,predictor)

All methods can be accessed in the same way; DA."test" where "test" is
the abbreviation given in the details of the testDA function.
