[![Build Status](https://travis-ci.org/hdsu-bioquant/ButchR.svg?branch=master)](https://travis-ci.org/hdsu-bioquant/ButchR)
[![codecov](https://codecov.io/gh/hdsu-bioquant/ButchR/branch/master/graph/badge.svg)](https://codecov.io/gh/hdsu-bioquant/ButchR)


<img src="vignettes/figs/ButchR/ButchR_logo.png" width="300">

Andres Quintero, Daniel Huebschmann & Sebastian Steinhauser  
15.07.2020  

`ButchR` is an R package providing functions to perform 
non-negative matrix factorization and postprocessing using TensorFlow.  

A detailed description 
of the software and an application to cells of the human hematopoietic system 
are available as a preprint: https://doi.org/10.1101/199547.  
Intermediate results for the analysis on hematopoietic cells are available on
zenodo: https://doi.org/10.5281/zenodo.800049.  


## How to install ButchR

### Install TensorFlow  

All the matrix decomposition algorithms implemented in ButchR run using 
TensorFlow > 2.0.0. Thus, a working installation of TensorFlow is needed to use 
the package.

There are several ways to install TensorFlow, for example using the R package 
`tensorflow`:

``` r
install.packages("tensorflow")
library(tensorflow)
install_tensorflow(version = "2.2.0")
```

Or, if there is a conda environment with TensorFlow installed, it can be 
activated before loading `ButchR`:

``` r
# It is important to set the environment before loading reticulate
reticulate::use_condaenv("tensorflow2env", required = TRUE)
library(reticulate)
py_config()
```

### Install ButchR   

ButchR can be installed by: 

``` r
remotes::install_github('wurst-theke/ButchR')
# same as devtools::install_github('wurst-theke/ButchR')
library(ButchR)
```

And a pre-build image to run ButchR inside RStudio can be pulled from Docker:
https://hub.docker.com/r/hdsu/butchr

`docker run --rm -p 8787:8787 -e USER=hdsu -e PASSWORD=pass hdsu/butchr`


## ShinyButchR

We also provide an interactive R/Shiny app `ShinyButchR` to perform NMF and 
explore the results interactively.
The live version can be used here: 
https://hdsu-bioquant.shinyapps.io/shinyButchR/

Or a pre-build image can be pulled from Docker:
https://hub.docker.com/r/hdsu/shinybutchr

`docker run --rm  -p 3838:3838 hdsu/shinybutchr`



## How to use ButchR


``` r
library(BiocStyle)
library(ButchR)
library(knitr)
library(ComplexHeatmap)
library(viridis)
library(tidyverse)
```

# Introduction


**NMF** (**nonnegative matrix factorization**) is a matrix decomposition
method. A description of the algorithm and it’s implementation can be
found e.g. in (Lee and Seung 1999). In 2003, Brunet et al. applied NMF
to gene expression data (Brunet et al. 2003). In 2010,
*[NMF](https://CRAN.R-project.org/package=NMF)*, an R package
implementing several NMF solvers was published (Gaujoux and Seoighe
2010). NMF basically solves the problem as illustrated in the following
figure (Image taken from
<a href="https://en.wikipedia.org/wiki/Non-negative_matrix_factorization" class="uri">https://en.wikipedia.org/wiki/Non-negative_matrix_factorization</a>):

![NMF](vignettes/figs/ButchR/NMF.png)

Here, *V* is an input matrix with dimensions *n* × *m*. It is decomposed
into two matrices *W* of dimension *n* × *l* and *H* of dimension
*l* × *m*, which when multiplied approximate the original matrix *V*.
*l* is a free parameter in NMF, it is called the factorization rank. If
we call the columns of *W* , then *l* corresponds to the number of
signatures. The decomposition thus leads to a reduction in complexity if
*l* \< *n*, i.e. if the number of signatures is smaller than the number
of features, as indicated in the above figure.

In 2015, Mejia-Roa et al. introduced an implementation of an NMF-solver
in CUDA, which lead to significant reduction of computation times by
making use of massive parallelisation on GPUs (Mejia-Roa et al. 2015).
Other implementations of NMF-solvers on GPUs exist.

It is the pupose of the package `ButchR` described here to provide
wrapper functions in R to these NMF-solvers in TensorFlow. Massive
parallelisation not only leads to faster algorithms, but also makes the
benefits of NMF accessible to much bigger matrices. Furthermore,
functions for estimation of the optimal factorization rank and post-hoc
feature selection are provided.

# The ButchR package


The matrix decomposition results are stored in an S4 object called
`ButchR_NMF`. `ButchR` provides functions to access the best
factorzation after *n* initailization *W* and *H* matrices for a given
factorzation rank.

A crucial step in data analysis with NMF is the determination of the
optimal factorization rank, i.e. the number of columns of the matrix *W*
or equivalently the number of rows of the matrix *H*. No consensus
method for an automatic evaluation of the optimal factorization rank has
been found to date. Instead, the decomposition is usually performed
iteratively over a range of possible factorization ranks and different
quality measures are computed for every tested factorization ranks. Many
quality measures have been proposed:

-   The `Frobenius reconstruction error`, i.e. the Frobenius norm of the
    residuals of the decomposition:
    \|\|*W* ⋅ *H* − *V*\|\|<sub>*F*</sub>

-   Criteria to assess the stability of the decomposition:

    -   The `cophenetic correlation coefficient`
    -   An `Amari type distance`
    -   `Silhouette values` over clusters of patterns extracted
        iteratively at the same factorization rank

The package `ButchR` provides a function to visualize all factorization
metrics.

# Example: leukemia data


Preparations

Load the example data

``` r
data(leukemia)
```

Now we are ready to start an NMF analysis.

NMF analysis
------------

### Call wrapper function

The wrapper function for the NMF solvers in the ButchR package is
`run_NMF_tensor`. It is called as follows:

``` r
k_min <- 2
k_max <- 4

leukemia_nmf_exp <- run_NMF_tensor(X = leukemia$matrix,
                                   ranks = k_min:k_max,
                                   method = "NMF",
                                   n_initializations = 10, 
                                   extract_features = TRUE)
```

    ## [1] "2020-07-16 17:50:42 CEST"
    ## Factorization rank:  2 
    ## [1] "NMF converged after  75,123,64,69,58,126,141,83,54,87 iterations"
    ## [1] "2020-07-16 17:50:42 CEST"
    ## Factorization rank:  3 
    ## [1] "NMF converged after  154,79,90,87,66,84,76,151,115,102 iterations"
    ## [1] "2020-07-16 17:50:44 CEST"
    ## Factorization rank:  4 
    ## [1] "NMF converged after  108,189,202,108,121,76,104,150,110,132 iterations"
    ## No optimal K could be determined from the Optimal K stat

Depending on the choice of parameters (dimensions of the input matrix,
number of iterations), this step may take some time. Note that the
algorithm updates the user about the progress in the iterations.

### normalize W matrix

To make the features in the *W* matrix comparable, the factorization is
normalized to make all columns of *W* sum 1.

``` r
leukemia_nmf_exp <- normalizeW(leukemia_nmf_exp)
```

Several functions to access the results are available:

### `HMatrix`

Returns the matrix `H` for the optimal decomposition (i.e. the one with
the minimal residual) for a specific factorization rank `k`. The number
of rows of the matrix `H` corresponds to the chosen factorization rank.

``` r
leukemia_Hk2 <- HMatrix(leukemia_nmf_exp, k = 2)
class(leukemia_Hk2)
```

    ## [1] "matrix" "array"

``` r
dim(leukemia_Hk2)
```

    ## [1]  2 38

``` r
kable(leukemia_Hk2[, 1:5])
```

|     ALL001|     ALL002|     ALL003|     ALL004|     ALL005|
|----------:|----------:|----------:|----------:|----------:|
|  1543455.8|  1356277.8|  1460730.1|  1420151.7|  1286558.2|
|   229239.7|   281053.9|   280561.9|   298372.2|   257627.2|

If no value for `k` is supplied, the function returns a list of
matrices, one for every factorization rank.

``` r
leukemia_Hlist <- HMatrix(leukemia_nmf_exp)
class(leukemia_Hlist)
```

    ## [1] "list"

``` r
length(leukemia_Hlist)
```

    ## [1] 3

``` r
kable(leukemia_Hlist$k2[, 1:5])
```

|     ALL001|     ALL002|     ALL003|     ALL004|     ALL005|
|----------:|----------:|----------:|----------:|----------:|
|  1543455.8|  1356277.8|  1460730.1|  1420151.7|  1286558.2|
|   229239.7|   281053.9|   280561.9|   298372.2|   257627.2|

### `WMatrix`

Returns the matrix `W` for the optimal decomposition (i.e. the one with
the minimal residual) for a specific factorization rank `k`. The number
of columns of the matrix `W` corresponds to the chosen factorization
rank.

``` r
leukemia_Wk2 <- WMatrix(leukemia_nmf_exp, k = 2)
class(leukemia_Wk2)
```

    ## [1] "matrix" "array"

``` r
dim(leukemia_Wk2)
```

    ## [1] 4452    2

``` r
kable(as.data.frame(leukemia_Wk2[1:5, ]))
```

|       |         V1|        V2|
|:------|----------:|---------:|
| A2M   |  0.0000453|  6.77e-05|
| AADAC |  0.0000793|  9.02e-05|
| AARS  |  0.0006768|  1.68e-04|
| ABAT  |  0.0001170|  1.41e-04|
| ABCA3 |  0.0000322|  8.74e-05|

If no value for `k` is supplied, the function returns a list of
matrices, one for every factorization rank.

``` r
leukemia_Wlist <- WMatrix(leukemia_nmf_exp)
class(leukemia_Wlist)
```

    ## [1] "list"

``` r
length(leukemia_Wlist)
```

    ## [1] 3

``` r
kable(as.data.frame(leukemia_Wlist$k2[1:5, ]))
```

|       |         V1|        V2|
|:------|----------:|---------:|
| A2M   |  0.0000453|  6.77e-05|
| AADAC |  0.0000793|  9.02e-05|
| AARS  |  0.0006768|  1.68e-04|
| ABAT  |  0.0001170|  1.41e-04|
| ABCA3 |  0.0000322|  8.74e-05|

### `FrobError`

Returns a data frame with as many columns as there are iterated
factorization ranks and as many rows as there are iterations per
factorization rank.

``` r
kable(FrobError(leukemia_nmf_exp))
```

|         k2|         k3|         k4|
|----------:|----------:|----------:|
|  0.5338607|  0.4778121|  0.4417868|
|  0.5342311|  0.4779813|  0.4408565|
|  0.5338451|  0.4783073|  0.4427480|
|  0.5338433|  0.4781617|  0.4410060|
|  0.5336356|  0.4780584|  0.4540046|
|  0.5335016|  0.4781451|  0.4413026|
|  0.5334857|  0.4779530|  0.4411095|
|  0.5335520|  0.4783140|  0.4411363|
|  0.5336552|  0.4779924|  0.4415452|
|  0.5335182|  0.4787099|  0.4409682|

Determine the optimal factorization rank
----------------------------------------

In NMF, Several methods have been described to assess the optimal
factorization rank. The ButchR package implements some of them:

-   *Frobenius error:* The most important information about the many
    iterated d ecompositions is the norm of the residual. In NMF this is
    often called the Frobenius error, as the Frobenius norm may be
    used.  
-   *Alexandrov Criterion:* In (Alexandrov et al. 2013) an approach is
    described in which a modified silhouette criterion is used to
    estimate the stability across iteration steps for one fixed
    factorization rank `k`.  
-   *Cophenetic correlation coefficient*  
-   *Amari distance*

The values of the computed factorization metrics can be accessed with
`OptKStats`:

``` r
kable(OptKStats(leukemia_nmf_exp))
```

<table style="width:100%;">
<colgroup>
<col style="width: 7%" />
<col style="width: 2%" />
<col style="width: 9%" />
<col style="width: 9%" />
<col style="width: 9%" />
<col style="width: 9%" />
<col style="width: 11%" />
<col style="width: 12%" />
<col style="width: 15%" />
<col style="width: 13%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">rank_id</th>
<th style="text-align: right;">k</th>
<th style="text-align: right;">min</th>
<th style="text-align: right;">mean</th>
<th style="text-align: right;">sd</th>
<th style="text-align: right;">cv</th>
<th style="text-align: right;">sumSilWidth</th>
<th style="text-align: right;">meanSilWidth</th>
<th style="text-align: right;">copheneticCoeff</th>
<th style="text-align: right;">meanAmariDist</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">k2</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">0.5334857</td>
<td style="text-align: right;">0.5337128</td>
<td style="text-align: right;">0.0002343</td>
<td style="text-align: right;">0.0004391</td>
<td style="text-align: right;">19.95919</td>
<td style="text-align: right;">0.9979595</td>
<td style="text-align: right;">0.9964766</td>
<td style="text-align: right;">0.0008694</td>
</tr>
<tr class="even">
<td style="text-align: left;">k3</td>
<td style="text-align: right;">3</td>
<td style="text-align: right;">0.4778121</td>
<td style="text-align: right;">0.4781435</td>
<td style="text-align: right;">0.0002538</td>
<td style="text-align: right;">0.0005307</td>
<td style="text-align: right;">29.94162</td>
<td style="text-align: right;">0.9980541</td>
<td style="text-align: right;">0.9765375</td>
<td style="text-align: right;">0.0008432</td>
</tr>
<tr class="odd">
<td style="text-align: left;">k4</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">0.4408565</td>
<td style="text-align: right;">0.4426464</td>
<td style="text-align: right;">0.0040295</td>
<td style="text-align: right;">0.0091031</td>
<td style="text-align: right;">37.69052</td>
<td style="text-align: right;">0.9422630</td>
<td style="text-align: right;">0.9591708</td>
<td style="text-align: right;">0.0237371</td>
</tr>
</tbody>
</table>

These quality measures can be displayed together:

### Generate plots to estimate optimal k

``` r
gg_plotKStats(leukemia_nmf_exp)
```

![](vignettes/figs/ButchR/kstats-1.png)

Visualize the matrix H (exposures)
----------------------------------

The matrices `H` may be visualized as heatmaps. We can define a meta
information object and annotate meta data:

``` r
heat_anno <- HeatmapAnnotation(df = leukemia$annotation[, c("ALL_AML", "Type")],
                               col = list(ALL_AML = c("ALL" = "grey80", 
                                                      "AML" = "grey20"),
                                          Type = c("-" = "white",
                                                   "B-cell" = "grey80",
                                                   "T-cell" = "grey20")))
```

And now display the matrices `H` with meta data annotation:

``` r
for(ki in k_min:k_max) {
  cat("\n")
  cat("  \n#### H matrix for k=",  ki, "   \n  ")
  #plot H matrix
  tmp_hmatrix <- HMatrix(leukemia_nmf_exp, k = ki)
  h_heatmap <- Heatmap(tmp_hmatrix,
                       col = viridis(100),
                       name = "Exposure",
                       clustering_distance_columns = 'pearson',
                       show_column_dend = TRUE,
                       top_annotation = heat_anno,
                       show_column_names = FALSE,
                       show_row_names = FALSE,
                       cluster_rows = FALSE)
  print(h_heatmap)
}
```

#### H matrix for k= 2

![](vignettes/figs/ButchR/hheatmap-1.png)

#### H matrix for k= 3

![](vignettes/figs/ButchR/hheatmap-2.png)

#### H matrix for k= 4

![](vignettes/figs/ButchR/hheatmap-3.png)

Feature selection
-----------------

### Row K-means to determine signature specific features

``` r
### Find representative regions.
# Get W for best K
leukemia_features <- SignatureSpecificFeatures(leukemia_nmf_exp,
                                               k = 4, 
                                               return_all_features = TRUE)
colnames(leukemia_features) <- paste0("Sign.", 1:4)
kable(head(leukemia_features))
```

|       |  Sign.1|  Sign.2|  Sign.3|  Sign.4|
|:------|-------:|-------:|-------:|-------:|
| A2M   |       1|       0|       0|       0|
| AADAC |       0|       0|       0|       1|
| AARS  |       0|       1|       1|       0|
| ABAT  |       0|       0|       0|       1|
| ABCA3 |       1|       0|       0|       1|
| ABCA4 |       0|       0|       0|       1|

### Feature visualization

``` r
# List of signature specific features
# leukemia_specific <- SignatureSpecificFeatures(leukemia_nmf_exp,
#                                                k = 4, 
#                                                return_all_features = FALSE)


leukemia_specific <- rownames(leukemia_features)[rowSums(leukemia_features) == 1]
leukemia_Wspecific <- WMatrix(leukemia_nmf_exp, k = 4)[leukemia_specific, ]
colnames(leukemia_Wspecific) <- paste0("Sign.", 1:4)

# normalize exposure score in W matrix across rows
leukemia_Wspecific <- leukemia_Wspecific/matrixStats::rowMaxs(leukemia_Wspecific)

# Display selected features on W matrix
w_heatmap <- Heatmap(leukemia_Wspecific,
                     col = inferno(100),
                     name = "W matrix",
                     clustering_distance_columns = 'pearson',
                     show_column_dend = TRUE,
                     show_column_names = TRUE,
                     show_row_names = FALSE,
                     cluster_rows = TRUE,
                     cluster_columns = FALSE)
w_heatmap

```

![](vignettes/figs/ButchR/wspecific-1.png)


# References


Alexandrov, LB, S Nik-Zainal, DC Wedge, SA Aparicio, S Behjati, AV
Biankin, GR Bignell, et al. 2013. “Signatures of Mutational Processes in
Cancer.” *Nature*. Nature Publishing Group.

Brunet, Jean-Philippe, Pablo Tamayo, Todd R. Golub, and Jill P. Mesirov.
2003. “Metagenes and Molecular Pattern Discovery Using Matrix
Factorization.” *PNAS*. PNAS.

Gaujoux, Renaud, and Cathal Seoighe. 2010. “A Flexible R Package for
Nonnegative Matrix Factorization.” *BMC Bioinformatics*. BMC.

Lee, Daniel D., and Sebastian Seung. 1999. “Learning the Parts of
Objects by Non-Negative Matrix Factorization.” *Nature*. Nature
Publishing Group.

Mejia-Roa, Edgardo, Daniel Tabas-Madrid, Javier Setoain, Carlos Garcia,
Francisco Tirado, and Alberto Pascual-Montano. 2015. “NMF-mGPU:
Non-Negative Matrix Factorization on Multi-GPU Systems.” *BMC
Bioinformatics*. BMC.
