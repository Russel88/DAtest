DAtest
======

This is a package for comparing different differential abundance methods
used in microbial marker-gene (e.g. 16S rRNA), RNA-seq and
protein/metabolite abundance analysis.

There are many methods for testing differential abundance and no gold
standard, but results can vary a lot between the different statistical
methods. The false positive rate and the power of the different methods
highly depends on the dataset. This package aims at aiding the analyst
in choosing a method for a specific dataset based on empirical testing.

#### The method goes as follows:

-   Shuffle predictor variable (E.g. case vs. control)
-   Spike in data for some randomly chosen features (OTU/gene/protein),
    such that they are associated with the shuffled predictor
-   Apply methods, and check:
    -   whether they can find the spike-ins
    -   whether the false positive rate is controlled

#### The workflow (details can be found below):

-   Compare methods with `testDA` function
    -   Check the results with `plot` or `summary`
    -   Choose method that has high AUC, and FPR not higher than ~0.05
-   Run data with the chosen test with DA."test" function, where "test"
    is the name of the test (see details with ?testDA)
    -   Check out your final results.

### Citation

Please cite the following publication if you use the DAtest package:

Russel *et al.* (2017) DAtest: A framework for comparing differential
abundance and expression methods. *Submitted*

Remember also to cite the method you end up using for your final
analysis (See [implemented methods](#implemented-methods) for links). If
there is no associated publication (e.g. for a t-test) you can cite as
follows: ."we used t-test as implemented in DAtest version x.x.x (Russel
*et al.* 2017)"

Overview of this tutorial
-------------------------

-   [Installation of packages](#installation-of-packages)
-   [How to compare methods](#how-to-compare-methods)
-   [How to run real data](#how-to-run-real-data)
-   [Implemented methods](#implemented-methods)
-   [Extra features](#extra-features)

Installation of packages
========================

It is advised not to have any packages loaded when installing DAtest.
Installation of DAtest and all dependencies has been tested to work on a
clean R version 3.4.1 by installing in the following order:

The DAtest package:

    install.packages("devtools")
    library(devtools)
    install_github("Russel88/DAtest")

#### The following are needed for *full* functionality

But the package will work without them

    source("https://bioconductor.org/biocLite.R")
    biocLite("DESeq2")
    biocLite("limma")
    biocLite("edgeR")
    biocLite("metagenomeSeq")
    biocLite("baySeq")
    biocLite("ALDEx2")
    biocLite("qvalue")

    biocLite("impute") # For SamSeq
    install.packages("samr") # SamSeq

    install.packages("pscl") # Zero-inflated Poisson and zero-inflated Negative Binomial
    install.packages("statmod") # For limma with blocking

-   RAIDA and ANCOM have to be installed from external sources:

<!-- -->

    # Dependencies (after the bioconductor packages are installed)
    install.packages(c("shiny","exactRankTests","openxlsx","protoclust","DT","coin","stringr"))

    # RAIDA:
    install.packages("https://cals.arizona.edu/~anling/software/RAIDA_1.0.tar.gz",repos = NULL)

    # ANCOM (by default not included, as it is very slow):
    download.file("https://www.niehs.nih.gov/research/resources/software/biostatistics/ancom/ancom_software.zip",destfile = "ANCOM.zip")
    unzip("ANCOM.zip",exdir=getwd())
    install.packages("ancom.R_1.1-3.tar.gz", repos = NULL)

**Note:** If installation fails for any of the bioconductor or external
packages (RAIDA, ANCOM) do not despair. `DAtest` will work seamlessly,
but will simply exclude methods that depends on these packages.

The following are suggested, but not needed:

    # For drawing Venn diagrams
    install.packages("eulerr")

    # For post-hoc testing (generalized) linear models
    install.packages("lsmeans")

How to compare methods
======================

Methods are compared with "False Positve Rate" (FPR), "Area Under the
(Receiver Operator) Curve" (AUC), and potentially also "Spike Detection
Rate" (Spike.detect.rate).

By shuffling the predictor variable and spiking a subset of the
features, we know a priori which features should be siginificant (the
spiked features) and which features shouldn't (the non-spiked features).

FPR indicates the proportion of non-spiked features that were found to
be significant (raw p-value below 0.05). With randomly generated data
the p-value distribution should be uniform (flat) and 5% of the raw
p-values are expected to be below 0.05. We therefore want an FPR at 0.05
or lower.

AUC is estimated by ranking all features by their respective p-value,
and finding the area under the ROC curve. For a good method p-values
from the spiked features should be among the lowest. For a perfect
method, with AUC = 1, all the lowest p-values are spiked features. An
AUC of 0.5 means that p-values from the spiked features are randomly
spread across the spectrum. An AUC below 0.5 means that p-values from
spiked features in general are *higher* than for non-spiked features.
Therefore, we want an AUC as high as possible; spiked features should
have low p-values.

Spike.detect.rate is the proportion of spiked features that are
significant after multiple correction of the p-values. It is therefore
the proportion of features you would expect to detect in a regular
analysis. The higher Spike.detect.rate, the better.

-   The intended workflow for choosing a method is:
    -   Omit methods with FPR significantly higher than 0.05
    -   Of remaining methods, choose the one with highest AUC
    -   If several methods have very similar AUCs, the Spike.detect.rate
        can be used to differentiate among the methods

### **Run the test:**

(if you have a phyloseq object, see details further down)

    mytest <- testDA(data, predictor, R = 10)

`data` is either a matrix or data.frame with taxa/genes/proteins as rows
and samples as columns (with rownames).

`predictor` is the predictor of interest, e.g. a factor denoting whether
samples are cases or controls (in the same order as columns in data).

`predictor` can be a factor with more than two levels, in which case
only the second level is spiked.

`predictor` can also be continuous/quantitative

**Important:** If `predictor` is numeric it will be treated as a
quantitative variable in the analyses. Therefore, if you, for example,
have two groups of samples and call the groups 0 and 1, it will be
treated as a quantitative variable and it will not run t-test,
wilcoxon-test, and others dedicated for two-class predictors. Write
`predictor` as `as.factor(predictor)` to treat it as a categorical
variable. The `testDA` function will produce a message telling you how
the `predictor` is read - Make sure this fits your expectations!

`R` denotes how many times the spike-in and FPR/AUC calculation should
be replicated. It is advised to use at least 10, but it can be set to 1
for a fast test of the function.

A subset of methods can be run by setting the `tests` argument.

#### **The function automatically uses multiple CPUs for fast execution**

The methods run in parallel, and by default the number of cores used is
one less than available. It has been tested to work on Windows, OS X and
Linux Debian. It can be changed with the `cores` argument. `cores = 1`
will turn off parallel computing. If the function is terminated before
ending you might get the following warning, which can safely be ignored:

    closing unused connection X (...)

But if you have terminated the function before it ended and your
computer runs slow, you might want to restart R to close all
connections.

If you run out of memory/RAM, reduce the number of cores used for
computing.

### *If you have more than 10k features (or 2k for paired analysis):*

Runtime of the different methods can vary quite a lot, and some methods
are simply unfeasible for datasets with several thousands of features.
The `runtimeDA` function can estimate the runtime of the different
methods on your dataset: It runs on small subsets of the data and then
extrapolates the runtime from these subsets. It will not be super
precise if you have tenths of thousands of features or more, but it
should give a decent estimate. The function prints how many minutes it
will take to run each method one time and `R` times. If you only use 1
core the actual runtime will be the sum of all the methods. With more
cores the runtime will decrease and approach the runtime of the slowest
method.

    runtimeDA(data, predictor)

### *If your predictor is categorical with more than two levels:*

There are generally two ways to output results with a categorical
predictor with multiple levels; either there is one p-value indicating
whether the categories are similar or different (e.g. in ANOVA), or
there is one p-value for each level, often where the first level is set
as intercept and the remaining levels are tested against the intercept
(e.g. in linear regression). For some methods you can choose which
option fits you, with other methods not, but it is crucial that this
option is similar in the `testDA` function as in the final analysis. Use
the `out.anova` argument to set this.

Below is a description of how the methods treat multi-class predictors:

All linear models (also GLMs) output results (including p-values) from
`anova`/`drop1` functions and are thus testing the `predictor` variable
in one go. For your final analysis you can run post-hoc tests for all
pairwise comparisons. If you are interested in testing treatments
against a common baseline/control (i.e. intercept), you can set
`out.anova = FALSE`. This will output results from the 2. level of the
`predictor` compared to the intercept. This is because only the 2. level
is spiked when `predictor` contains multiple levels. In your final
analysis you can get an output with all p-values (See more
[here](#how-to-run-real-data)).

All limma models output results (including p-values) from `topTable`
testing all levels (minus 1) against the intercept. This can be changed
with `out.anova`. `out.anova = FALSE` will output results from the 2.
level of the predictor. In your final analysis you can set
`allResults = TRUE` and use `topTable` on the output to get the desired
results.

DESeq2 is set to run Log Ratio Test (LRT) and is thus testing all levels
of the `predictor` in one go.

EdgeR is set to test if all levels of `predictor` (minus intercept) are
zero.

### *If you have a paired/blocked experimental design:*

E.g. repeated samples from same patients, or control/treatments inside
blocks.

The `paired` argument should be a factor with the ID of the
patient/subject/block (in the same order as columns in `data`):

    mytest <- testDA(data, predictor, paired = subjectID)

When a `paired` argument is provided, the `predictor` is shuffled within
the levels of the `paired` factor.

Paired analysis can be very slow. If you simply can't wait to see the
results remove "neb" from the `tests` argument.

**Some details on how methods use the paired argument:**

t-test, wilcox test, friedman test, quade test, SAMseq and permutation
test expect a balanced "unreplicated" design with one value for each
combination of levels in the `paired` and `predictor` variables.

Negbinom/poisson glm, linear regression and limma models use the
`paired` variable as a random intercept (i.e. they become mixed-effect
models) and are very flexible regarding design and work well for
unbalanced designs.

EdgeR, DESeq2 and metagenomeSeq ZIG include the `paired` variable as a
covariable in the model and are also generally flexible in the design.

ANCOM use the `paired` variable in a repeated measures manner

### *If you have non-relative abundances, e.g. for normalized protein abundance:*

    mytest <- testDA(data, predictor, relative = FALSE)

If `relative=FALSE` the values in `data` will NOT be normalized by the
sum of the samples for "ttt", "ltt", "ltt2", "wil", "per", "lrm", "llm",
"llm2", "lim", "lli", "lli2", "pea", "spe", "aov", "lao", "lao2" and
"kru", and there will NOT be an offset of log(librarySize) for "neb",
"poi", "qpo", "zpo" and "znb". Furthermore, tests that are specifically
designed to handle sequence data are turned off: DESeq2, EdgeR,
MetagenomeSeq, ANCOM, RAIDA, ALDEx2 and baySeq. Therefore, the `data` is
analysed as provided, except for the tests that log-transform the data:
"ltt", "llm", "lli" and "lao".

### *If you have covariates:*

The `covars` argument should be a named list with the covariables (in
the same order as columns in `data`):

    mytest <- testDA(data, predictor, covars = list(Age = subject.age, Date = exp.date))

### *If you have a phyloseq object:*

`data` can also be a phyloseq object. In this case, the `predictor`,
`paired` and `covars` arguments are the names of the variables in
`sample_data(data)`:

    mytest <- testDA(data, predictor = "Time", paired = "Patient", covars = c("Age","ExpDate"))

**Check the output:**
---------------------

    plot(mytest)
    summary(mytest)

    # Details from the run:
    mytest$details

    # Average run times of each method:
    mytest$run.times

**Note:**

As ANCOM and SAMseq do not output p-values, AUC for "anc" and "sam" is
based on pseudo p-values. They are calculated from the statistics/scores
as these are perfectly ranked according to detection/significance
calling. FPR is also based on pseudo p-values for "anc" and "sam", but
as these cannot be adjusted as nominal p-values, FPR for these methods
is the final false discovery rate and we should expect an FPR close to 0
for these two methods, unless you are willing to accept some false
positives. This can be tuned with the `sig`/`multcorr` ("anc") and
`fdr.output` ("sam") arguments.

P-values for baySeq are defined as 1 - posterior likelihoods.

### *How to choose the effect size:*

The default effect size of 2 (equal to a log2 fold change of 1) might
not fit for your data. If the best method has an AUC below 0.7 you might
want to increase the effect size. In contrast, if several methods have
an AUC close to 1, you might want to decrease the effect size (towards
one) to better differentiate the methods. It is also possible to test
for negative associations with effect sizes between 0 and 1, but this
has not been tested extensively.

How to run real data
====================

All tests can easily be run with the original data. E.g. edgeR exact
test:

    res.ere <- DA.ere(data, predictor)

All methods can be accessed in the same way; DA."test" where "test" is
the abbreviation given in the details of the `testDA` function.

It is advised to set `allResults = TRUE` for checking final results. For
all methods where relevant, this will output the raw results, often in a
list with each element corresponding to a feature (OTU/gene/protein).
For published methods, it is advised to check their tutorials on how to
read to output.

-   ***IMPORTANT:***
    -   Set `out.anova` similar in `testDA` as in your final analysis
        for reliable comparison
    -   If your `predictor` has more than two levels you might have to
        set the `by` argument for "zig" (this is by default = 2)
    -   All linear models (including all GLMs) output the p-value from
        an `anova`/`drop1` function. This can be changed with the
        `out.anova` argument
    -   All limma models (lim,lli,lli2,vli) test all levels of the
        `predictor` against an intercept (with `topTable`). This can be
        changed with the `out.anova` argument
    -   For ANCOM: If the FPR = 0, you would not expect false positives
        with the default settings. If "anc" has an FPR &gt; 0, set
        `multcorr = 1`

For linear models the `drop1`/`anova` functions can be used to test
significance of the `predictor` and `covars` variables:

    ### Apply drop1 for each feature and output the adjusted p-values:
    # Works on "zpo", "znb", "qpo", "neb", "poi". Non-paired "lrm", "llm", "llm2"
    results <- DA.lrm(data, predictor, allResults = TRUE)
    res.drop1 <- DA.drop1(results)

    ### Apply anova for each feature and output the adjusted p-values:
    # Works on "lrm", "llm", "llm2". Non-paired "neb"
    results <- DA.lrm(data, predictor, allResults = TRUE)
    res.anova <- DA.anova(results)

For anova and linear models we can also run post-hoc tests for all
pairwise comparisons of multi-class `predictor`/`covars`.

    ### Apply TukeyHSD for each feature for a selected variable and output the adjusted p-values:
    # Works on "aov", "lao", "lao2"
    results <- DA.aov(data, predictor, allResults = TRUE)
    res.tukey <- DA.TukeyHSD(results, variable = "predictor") # variable can also be the name of a covar

    ### Apply lsmeans for each feature for a selected variable and output the adjusted p-values:
    # This requires the lsmeans package.
    # Works on "poi", "neb", "lrm", "llm", "llm2", "qpo", "znb", "zpo"
    results <- DA.lrm(data, predictor, allResults = TRUE)
    res.lsm <- DA.lsmeans(results, variable = "predictor") # variable can also be the name of a covar

    # For paired "lrm", "llm", "llm2" the original predictor variable has to be supplied. For example:
    results <- DA.lrm(data, predictor = mypred, paired = SubjectID,  allResults = TRUE)
    res.lsm <- DA.lsmeans(results, variable = "predictor", predictor = mypred)

    # and if covars are used, they also need to be supplied for paired "lrm", "llm", "llm2". For example:
    results <- DA.lrm(data, predictor = mypred, paired = SubjectID, covars = list(covar1 = mycovar),  allResults = TRUE)
    res.lsm <- DA.lsmeans(results, variable = "predictor", predictor = mypred, covars = list(covar1 = mycovar))

**Alternatively, run all (or several) methods and check which features
are found by several methods.**

    # Run many methods:
    res.all <- allDA(data, predictor)

    # Adjusted p-values from all methods (detection/no-detection from sam and anc)
    res.all$adj

    # Estimates/fold.changes from all methods which output anything relevant for this
    res.all$est

    # Venn (Euler) diagram of detected features from selected methods:
    # This requires the eulerr package
    vennDA(res.all, tests = c("wil","ttt","ltt"))

    # Split Venn (Euler) diagram in significant features with either positive or negative fold changes (if possible)
    vennDA(res.all, tests = c("wil","ttt","ltt"), split = TRUE)

    # See results from a method (e.g. t.test "ttt"):
    View(res.all$results$ttt)

A subset of methods can be run by setting the `tests` argument. E.g.
only those performing well based on results from `testDA`.

Implemented methods
===================

### Is your favorite method missing?

Either add it yourself [(see under 'Extra features')](#extra-features),
or write to me, preferably with a code snippet of the implementation
(see email in Description).

-   per - [Permutation test with user defined test
    statistic](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8)
    (See below for description of the paired permutation test)
-   bay -
    [baySeq](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-422)
-   adx - [ALDEx t-test and
    wilcoxon](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0067019)
-   wil - [Wilcoxon Rank
    Sum](https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test) on
    relative abundances (Paired is a [Wilcoxon Signed Rank
    test](https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test)
-   ttt - [Welch t.test](https://en.wikipedia.org/wiki/Welch%27s_t-test)
    on relative abundances (Paired is [paired
    t-test](https://en.wikipedia.org/wiki/Student%27s_t-test#Dependent_t-test_for_paired_samples))
-   ltt - [Welch
    t.test](https://en.wikipedia.org/wiki/Welch%27s_t-test), but reads
    are first transformed with log(abundance + delta) then turned into
    relative abundances
-   ltt2 - [Welch
    t.test](https://en.wikipedia.org/wiki/Welch%27s_t-test), but with
    relative abundances transformed with log(relative abundance + delta)
-   neb - [Negative binomial
    GLM](https://en.wikipedia.org/wiki/Negative_binomial_distribution)
    (The paired version is a [mixed-effect
    model](https://en.wikipedia.org/wiki/Mixed_model))
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
    (The paired version is a model with the paired variable
    as covariate)
-   ds2 -
    [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)
    (The paired version is a model with the paired variable
    as covariate)
-   lim -
    [LIMMA](https://link.springer.com/chapter/10.1007%2F0-387-29362-0_23?LI=true)
    (The paired version is using the block argument in lmFit)
-   lli -
    [LIMMA](https://link.springer.com/chapter/10.1007%2F0-387-29362-0_23?LI=true),
    but reads are first transformed with log(abundance + delta) then
    turned into relative abundances
-   lli2 -
    [LIMMA](https://link.springer.com/chapter/10.1007%2F0-387-29362-0_23?LI=true),
    but with relative abundances transformed with log(relative
    abundance + delta)
-   vli - [LIMMA with
    voom](https://link.springer.com/chapter/10.1007%2F0-387-29362-0_23?LI=true)
    (The paired version is using the block argument in lmFit)
-   kru - [Kruskal-Wallis
    test](https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance)
    on relative abundances
-   aov - [ANOVA](https://en.wikipedia.org/wiki/Analysis_of_variance) on
    relative abundances
-   lao - [ANOVA](https://en.wikipedia.org/wiki/Analysis_of_variance),
    but reads are first transformed with log(abundance + delta) then
    turned into relative abundances
-   lao2 - [ANOVA](https://en.wikipedia.org/wiki/Analysis_of_variance),
    but with relative abundances transformed with log(relative
    abundance + delta)
-   lrm - [Linear
    regression](https://en.wikipedia.org/wiki/Linear_regression) on
    relative abundances (The paired version is a [mixed-effect
    model](https://en.wikipedia.org/wiki/Mixed_model))
-   llm - [Linear
    regression](https://en.wikipedia.org/wiki/Linear_regression), but
    reads are first transformed with log(abundance + delta) then turned
    into relative abundances
-   llm2 - [Linear
    regression](https://en.wikipedia.org/wiki/Linear_regression), but
    with relative abundances transformed with log(relative abundance +
    delta)
-   rai -
    [RAIDA](https://academic.oup.com/bioinformatics/article/31/14/2269/256302/A-robust-approach-for-identifying-differentially?searchresult=1)
-   spe - [Spearman Rank
    Correlation](https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient)
-   pea - [Pearson
    Correlation](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient)
-   poi - [Poisson
    GLM](https://en.wikipedia.org/wiki/Poisson_distribution) (The paired
    version is a [mixed-effect
    model](https://en.wikipedia.org/wiki/Mixed_model))
-   qpo - [Quasi-poisson
    GLM](https://en.wikipedia.org/wiki/Quasi-likelihood)
-   zpo - [Zero-inflated Poisson
    GLM](https://en.wikipedia.org/wiki/Zero-inflated_model)
-   znb - [Zero-inflated Negative Binomial
    GLM](https://en.wikipedia.org/wiki/Zero-inflated_model)
-   fri - [Friedman Rank Sum
    test](https://en.wikipedia.org/wiki/Friedman_test)
-   qua - [Quade test](http://rcompanion.org/handbook/F_11.html)
-   anc - [ANCOM](https://www.ncbi.nlm.nih.gov/pubmed/26028277) (by
    default not included, as it is very slow)
-   sam - [SAMseq](http://statweb.stanford.edu/~tibs/SAM/)

### Paired permutation test

A paired permutation test is implemented specifically for this package.
The test is similar to the original, but with a different test statistic
and permutation scheme. The permutations are constrained in the paired
version such that the predictor is only permuted within each level of
the paired argument (e.g. subjects). The test statistic first finds the
log-ratio between the two predictor levels (e.g. case and control) for
each level of the paired argument and the final statistic is the mean of
these log-ratios.

Extra features
==============

#### Adding user-defined methods

"zzz" (`DA.zzz`) is a placeholder for user-defined methods. You have to
supply it with a function whose input is: A `count_table` (data.frame,
samples are columns, rownames indicate features), a `predictor`
(vector), a `paired` variable (factor), and a `covars` argument (named
list with vectors). It is OK if your function doesn't use
`paired`/`covars`, but they have to be there in the arguments. The
output from the user-defined function should be a data.frame that at
least includes: The names of the features ("Feature"), p-values
("pval"), and name of method ("Method").

See example below on how to include a simple t-test on relative
abundances:

    # Define our function
    myfun <- function(count_table, predictor, paired, covars){ # These fours arguments should not be altered
      
      # Relative abundance
      rel <- apply(count_table, 2, function(x) x/sum(x))
      
      # t-test function
      tfun <- function(x){
        tryCatch(t.test(x ~ predictor)$p.value, error = function(e){NA}) 
      }
      
      # P-values for each feature
      pvals <- apply(rel, 1, tfun)
      
      # Collect and return data
      df <- data.frame(Feature = rownames(count_table),
                       pval = pvals)
      df$Method <- "My own t-test"
      return(df)
      
    }

    # Test that it works on real data
    res <- DA.zzz(data, predictor, FUN = myfun)

    # Use it in the testDA function
    mytest <- testDA(data, predictor, tests = c("zzz","ttt","wil"), args = list(zzz = list(FUN = myfun)))

    # Several user-defined methods can be included by simply putting numbers after "zzz" in both the "tests" and "args" arguments:
    mytest <- testDA(data, predictor, tests = c("zzz","zzz2","ttt","wil"), args = list(zzz = list(FUN = myfun),
                                                                                       zzz2 = list(FUN = myfun2)))

**Caution:** If "zzz" is in the `tests` argument, there will be no
checks of whether any of the supplied methods are suitable for the data.
If e.g. your `predictor` is quantitative and "ttt" is in the `tests`
argument, it will try to run a t-test and fail miserably.

#### Results from all the runs:

    print(mytest)

#### See the output from the individual methods. E.g. "ere" first run:

    View(mytest$results[[1]]["ere"])

#### Plot p-value histograms for all the non-spiked features:

    plot(mytest, p = TRUE)

### Passing arguments to the different tests

Additional arguments can be passed to the internal functions with the
`args` argument. It should be structured as a list with elements named
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
-   ere - Passed to calcNormFactors, estimateCommonDisp,
    estimateTagwiseDisp and exactTest
-   erq - Passed to calcNormFactors, estimateDisp, glmQLFit and
    glmQLFTest
-   msf - Passed to fitFeatureModel
-   zig - Passed to fitZig
-   ds2 - Passed to DESeq
-   lim - Passed to eBayes and lmFit
-   lli - Passed to eBayes, lmFit and DA.lli
-   lli2 - Passed to eBayes, lmFit and DA.lli2
-   kru - Passed to kruskal.test
-   aov - Passed to aov
-   lao - Passed to aov and DA.lao
-   lao2 - Passed to aov and DA.lao2
-   lrm - Passed to lm and lme
-   llm - Passed to lm, lme and DA.llm
-   llm2 - Passed to lm, lme and DA.llm2
-   rai - Passed to raida
-   spe - Passed to cor.test
-   pea - Passed to cor.test
-   poi - Passed to glm and glmer
-   qpo - Passed to glm
-   vli - Passed to voom, eBayes and lmFit
-   zpo - Passed to zeroinfl
-   znb - Passed to zeroinfl
-   fri - Passed to friedman.test
-   qua - Passed to quade.test
-   anc - Passed to ANCOM
-   sam - Passed to SAMseq
