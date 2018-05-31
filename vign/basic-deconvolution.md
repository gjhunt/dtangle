A Basic Deconvolution Example
-----------------------------

In this vignette we will work through a simple example of deconvolving cell type proportions from DNA microarray data. We work with a data set created from rats and introduced by [Shen-Orr et al](https://www.nature.com/nmeth/journal/v7/n4/abs/nmeth.1439.html). This is available on GEO with accession [GSE19830](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19830). The data set we will work with is subset of the Shen-Orr data and is included in the `dtangle` package under the name `shen_orr_ex`. Alternatively, we can access this and other data sets data set through the supplementary `dtangle.data` package we have made available [here](https://umich.box.com/v/dtangledatapkg). More information about the data set is available as part of the `R` help, `?shen_orr_ex`. First load up the data set.

``` r
library('dtangle')
names(shen_orr_ex)
```

    ## [1] "data"       "annotation" "name"

In this data set rat brain, liver and lung cells have been mixed together in various proportions the resulting mixtures were analyzed with DNA microarrays. The mixing proportions are encoded in the mixture matrix

``` r
truth = shen_orr_ex$annotation$mixture
head(truth)
```

    ##           Liver Brain Lung
    ## GSM495209     1     0    0
    ## GSM495210     1     0    0
    ## GSM495211     1     0    0
    ## GSM495212     0     1    0
    ## GSM495213     0     1    0
    ## GSM495214     0     1    0

Each row of this matrix is a sample and each column gives the mixing proportions of the cell types in each sample. From this we can extract out the pure samples of each of the three cell types.

``` r
pure_samples <- lapply(1:3, function(i) {
    which(truth[, i] == 1)
})
names(pure_samples) = colnames(truth)
pure_samples
```

    ## $Liver
    ## GSM495209 GSM495210 GSM495211 
    ##         1         2         3 
    ## 
    ## $Brain
    ## GSM495212 GSM495213 GSM495214 
    ##         4         5         6 
    ## 
    ## $Lung
    ## GSM495215 GSM495216 GSM495217 
    ##         7         8         9

The RMA-summarized gene expression data generated as part of the Shen-Orr experiment is accessible under `data$log`,

``` r
Y <- shen_orr_ex$data$log
Y[1:4,1:4]
```

    ##           X1367452_at X1367453_at X1367454_at X1367455_at
    ## GSM495209    8.869687    8.757657    9.320518   10.070805
    ## GSM495210    8.938069    8.741337    9.225232   10.096493
    ## GSM495211    8.883494    8.714461    9.302988   10.089371
    ## GSM495212    9.918893    8.974574    8.293384    9.363512

Each row is a different individual and each column is a particular gene. The values of the matrix are log<sub>2</sub> RMA-summarized gene expressions.

The first step in running `dtangle` is to identify marker genes for each cell type. These may be provided by the scientist if they are already known or may be determined by `dtangle` or another algorithm. To find marker genes using `dtangle` we pass the following arugments to the `find_markers` function:

1.  the data matrix, `Y`,

2.  the list of pure samples for each type, `pure_samples`,

3.  the data type, `data_type`,

4.  the method used to rank markers, `marker_method`.

``` r
marker_list = find_markers(Y,pure_samples,data_type="microarray-gene",marker_method='ratio')
```

The function `find_markers` will determine which genes to use as markers of each cell type. The function returns a list of two elements. The first element is `L` which is a list where the *i*<sup>*t**h*</sup> element is a vector of marker indices (columns of *Y*) for the *i*<sup>*t**h*</sup> type ranked in decreasing order of utility (best markers listed first) according to the chosen method.

``` r
lapply(marker_list$L,head)
```

    ## [[1]]
    ## [1] 466 597 353 196 387 537
    ## 
    ## [[2]]
    ## [1] 400 593 394 541 479 583
    ## 
    ## [[3]]
    ## [1] 115 210 163 213 395 117

The second element of the list returned by `find_markers` is a list `V` (with the same structure as `L`). For each of the indices in `L` the list `V` contains the marker score as determined by the method utilized in `find_markers`. Larger values are better. The meaning of this value in `V` depends on the choice of `marker_method`. For `marker_method=ratio` (the default method) the values in `V` are the ratio of the estimated amount the particular gene is expressed in each cell type over the sum of the expression of this gene in all other cell types.

``` r
lapply(marker_list$V,head)
```

    ## [[1]]
    ## [1] 142.27505 128.57801 109.90386  49.51434  45.53533  33.57015
    ## 
    ## [[2]]
    ## [1] 211.60346 105.44800  99.10636  98.68013  86.97392  71.04771
    ## 
    ## [[3]]
    ## [1] 599.35967  54.71688  21.36921  18.41595  17.92601  12.68834

After we have ranked our markers with `find_markers` we need to determine how many markers to use for each cell type. The simplest way to do this is to choose, say, the top 10% of all marker genes for each type.

``` r
q = .1
quantiles = lapply(marker_list$V,function(x)quantile(x,1-q))
K = length(pure_samples)
n_choose = sapply(1:K,function(i){max(which(marker_list$V[[i]] > quantiles[[i]]))})
n_choose
```

    ## [1] 23 17 22

Now that we have ranked the genes as markers for each type and chosen how many marker genes to use for each cell type we can run the `dtangle` deconvolution algorithm.

``` r
marks = marker_list$L
dc <- dtangle(Y,pure_samples,n_choose,data_type='microarray-gene',markers=marks)
```

providing to the `dtangle` function the arguments:

1.  `Y`, our microarray data

2.  the list of `pure_samples`

3.  the number of markers to use for each cell type, `n_choose`

4.  the `data_type`

5.  the list of ranked markers (output from `find_markers`) to the `markers` argument.

The `dtangle` function returns to us a list with elements

1.  `estimates`, the estimated mixing proportions for each type for each sample

2.  `markers`, the markers used for each cell type

3.  `n_choose`, how many markers we used for each cell type

4.  `gamma`, the value of the sensitivity parameter used.

``` r
dc
```

    ## $estimates
    ##                  [,1]        [,2]       [,3]
    ## GSM495209 0.986313938 0.002472217 0.01121384
    ## GSM495210 0.984989195 0.002588476 0.01242233
    ## GSM495211 0.985869822 0.002804773 0.01132541
    ## GSM495212 0.004903505 0.984718617 0.01037788
    ## GSM495213 0.004815152 0.984502313 0.01068254
    ## GSM495214 0.004985250 0.983932945 0.01108181
    ## GSM495215 0.005576067 0.002590890 0.99183304
    ## GSM495216 0.005421071 0.002734625 0.99184430
    ## GSM495217 0.005672315 0.002836182 0.99149150
    ## GSM495218 0.046902166 0.241985915 0.71111192
    ## GSM495219 0.050212583 0.234180104 0.71560731
    ## GSM495220 0.049669746 0.240352467 0.70997779
    ## GSM495221 0.701046642 0.047015320 0.25193804
    ## GSM495222 0.689919421 0.047778621 0.26230196
    ## GSM495223 0.696489842 0.043637242 0.25987292
    ## GSM495224 0.191040239 0.769652721 0.03930704
    ## GSM495225 0.195784096 0.764120085 0.04009582
    ## GSM495226 0.195181583 0.765951033 0.03886738
    ## GSM495227 0.658921776 0.293530022 0.04754820
    ## GSM495228 0.664760588 0.286800243 0.04843917
    ## GSM495229 0.660824210 0.290859187 0.04831660
    ## GSM495230 0.403629834 0.522438531 0.07393163
    ## GSM495231 0.405785448 0.517112943 0.07710161
    ## GSM495232 0.413723321 0.509975731 0.07630095
    ## GSM495233 0.547065375 0.210951907 0.24198272
    ## GSM495234 0.548696486 0.213566737 0.23773678
    ## GSM495235 0.547857533 0.209863933 0.24227853
    ## GSM495236 0.482641393 0.344216545 0.17314206
    ## GSM495237 0.486649434 0.343037127 0.17031344
    ## GSM495238 0.491947779 0.337544956 0.17050726
    ## GSM495239 0.534081406 0.338658027 0.12726057
    ## GSM495240 0.532174471 0.343907982 0.12391755
    ## GSM495241 0.529859974 0.342166427 0.12797360
    ## GSM495242 0.470093476 0.451703268 0.07820326
    ## GSM495243 0.464890547 0.456273497 0.07883596
    ## GSM495244 0.470559785 0.449350598 0.08008962
    ## GSM495245 0.553621829 0.401976274 0.04440190
    ## GSM495246 0.552694400 0.403056873 0.04424873
    ## GSM495247 0.558667273 0.394950337 0.04638239
    ## GSM495248 0.596594442 0.384387861 0.01901770
    ## GSM495249 0.596208829 0.384415005 0.01937617
    ## GSM495250 0.592328501 0.388616006 0.01905549
    ## 
    ## $markers
    ## $markers[[1]]
    ##  [1] 466 597 353 196 387 537 420 570 434 436 307 543 105 303 324 549 486
    ## [18] 458 595 347 565 488 594
    ## 
    ## $markers[[2]]
    ##  [1] 400 593 394 541 479 583 526 467 498 384 505 519 348 502 311 413 361
    ## 
    ## $markers[[3]]
    ##  [1] 115 210 163 213 395 117 119 365 133 250 214 126 448 228 429 123 399
    ## [18] 461 415 489 501 334
    ## 
    ## 
    ## $n_choose
    ## [1] 23 17 22
    ## 
    ## $gamma
    ## [1] 0.6999978

We can plot our estimates against the known truth as follows

``` r
phats <- dc$estimates
plot(truth,phats,xlab="Truth",ylab="Estimates",xlim=c(0,1),ylim=c(0,1))
abline(coef=c(0,1))
```

![](basic-deconvolution_files/figure-markdown_github/unnamed-chunk-11-1.png)

If desired, we can specify the value of the sensivity parameter `gamma` numerically instead of letting `dtangle` choose it based upon the `data_type`. For example,

``` r
dc <- dtangle(Y,pure_samples,n_choose,gamma=.7,markers=marks)
phats <- dc$estimates
plot(truth,phats,xlab="Truth",ylab="Estimates",xlim=c(0,1),ylim=c(0,1))
abline(coef=c(0,1))
```

![](basic-deconvolution_files/figure-markdown_github/unnamed-chunk-12-1.png)

We can view the pre-selected values for `gamma` for each data type by using the function `get_gamma`

``` r
get_gamma('microarray-probe')
```

    ## [1] 0.4522564

``` r
get_gamma('microarray-gene')
```

    ## [1] 0.6999978

``` r
get_gamma('rna-seq')
```

    ## [1] 0.9433902

We can also specify the number of markers to be the same for each cell type by providing a single number to `n_choose`

``` r
dc <- dtangle(Y,pure_samples,n_choose=5,data_type='microarray-gene',markers=marks)
phats <- dc$estimates
plot(truth,phats,xlab="Truth",ylab="Estimates",xlim=c(0,1),ylim=c(0,1))
abline(coef=c(0,1))
```

![](basic-deconvolution_files/figure-markdown_github/unnamed-chunk-14-1.png)

Alternatively, we can manually specify the number of markers to use for each cell type manually,

``` r
dc <- dtangle(Y,pure_samples,n_choose=c(5,6,7),data_type='microarray-gene',markers=marks)
phats <- dc$estimates
plot(truth,phats,xlab="Truth",ylab="Estimates",xlim=c(0,1),ylim=c(0,1))
abline(coef=c(0,1))
```

![](basic-deconvolution_files/figure-markdown_github/unnamed-chunk-15-1.png)

We can test different methods of choosing markers by specifying the `marker_method` argument. Notice that if we don't calculate the markers in advance (e.g. by using `find_markers`) then `dtangle` handles the markers internally by using `find_markers` with `method='ratio'`. A description of the marker choosing methods can be found in the help pages for `dtangle`.

``` r
dc <- dtangle(Y, pure_samples,n_choose,data_type='microarray-gene',marker_method = 'ratio')
phats <- dc$estimate
plot(truth,phats,xlab="Truth",ylab="Estimates",xlim=c(0,1),ylim=c(0,1))
abline(coef=c(0,1))

dc2 <- dtangle(Y, pure_samples,n_choose,data_type='microarray-gene',marker_method = 'diff')
phats2 <- dc2$estimates
points(truth,phats2,col='blue')

dc3 <- dtangle(Y, pure_samples,n_choose,data_type='microarray-gene',marker_method = 'regression')
phats3 <- dc3$estimates
points(truth,phats3,col='red')
```

![](basic-deconvolution_files/figure-markdown_github/unnamed-chunk-16-1.png)
