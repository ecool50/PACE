# Comparing proximity effects between conditions

## Introduction

The *[PACE](https://bioconductor.org/packages/3.23/PACE)* manual fits a
single section and asks how expression changes with proximity to each
neighbouring cell type. This vignette asks the next question: **does
that proximity effect differ between two groups of patients?**

Answering it needs a cohort rather than a section, and it adds one
component to the model. Alongside the spatial cell state block, which is
the proximity effect common to everyone, PACE estimates a **responder
spatial state** block: the part of the proximity effect that differs
between conditions. Read the two together, a gene can respond to a
neighbour in both arms equally (spatial, no responder signal), or
respond differently depending on outcome (responder signal).

## The melanoma subset

The package ships a subset of a CosMx melanoma cohort of patients
profiled after immune checkpoint blockade, contrasting progressive
disease (PD) against non-progressive disease (non-PD, pooling stable,
partial and complete responders).

``` r

library(PACE)
library(SpatialExperiment)
library(ggplot2)
library(dplyr)

spe <- readRDS(system.file("extdata", "mel_cosmx_subset.rds", package = "PACE"))
spe
#> class: SpatialExperiment 
#> dim: 927 8454 
#> metadata(0):
#> assays(1): counts
#> rownames(927): AATK ABL1 ... YES1 ZFP36
#> rowData names(0):
#> colnames(8454): Cell53068 Cell53069 ... Cell21330 Cell21344
#> colData names(4): cellType Responder image sample_id
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> spatialCoords names(2) : x y
#> imgData names(0):
table(spe$cellType)
#> 
#>      B_Cell Endothelial  Fibroblast  Macrophage      T_Cell      Tumour 
#>         174         182         211         453         437        6997
table(unique(data.frame(image = spe$image, arm = spe$Responder))$arm)
#> 
#> nonPD    PD 
#>    19     7
```

One section per patient, so the condition contrast is between patients,
not within a section. Each section here is a contiguous spatial crop of
the published one; the derivation is in
`inst/scripts/make-mel-cosmx-subset.R`.

> **Provenance.** These data are a modified subset of the deposit
> accompanying Dong, Su, Kluger, Fan and Kluger, *SIMVI disentangles
> intrinsic and spatial-induced cellular states in spatial omics data*
> ([doi:10.5281/zenodo.14708000](https://doi.org/10.5281/zenodo.14708000)),
> used here under CC-BY-4.0.

> **Read this before drawing conclusions.** The subset is 8,454 of the
> cohort’s 56,274 cells, sized so the package stays small enough to
> distribute. It is an illustration of the workflow, not an analysis.
> The section on statistical power below shows what that costs, with
> numbers.

## Fitting with a condition

The call is the single-section one plus three arguments:

- `condition_col` names the column holding the contrast, and switches on
  the responder spatial state block.
- `kernel_per_image = TRUE` builds neighbourhoods within each section,
  so no cell is ever a neighbour of a cell from another patient.
- `image_re = "intercept"` gives each section a random intercept, so
  patient-level differences in overall expression are absorbed rather
  than being attributed to the condition.

``` r

fit <- paceFit(spe,
               celltype_col     = "cellType",
               condition_col    = "Responder",
               image_col        = "image",
               kernel_per_image = TRUE,
               image_re         = "intercept",
               dispersion       = "nb1",
               verbose          = FALSE)
#>  - Computing 927 x 78 likelihood matrix.
#>  - Likelihood calculations took 0.06 seconds.
#>  - Fitting model with 78 mixture components.
#>  - Model fitting took 0.05 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.01 seconds.
#>  - Computing 927 x 254 likelihood matrix.
#>  - Likelihood calculations took 0.18 seconds.
#>  - Fitting model with 254 mixture components.
#>  - Model fitting took 0.12 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>  - Computing 927 x 265 likelihood matrix.
#>  - Likelihood calculations took 0.19 seconds.
#>  - Fitting model with 265 mixture components.
#>  - Model fitting took 0.11 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>  - Computing 927 x 78 likelihood matrix.
#>  - Likelihood calculations took 0.06 seconds.
#>  - Fitting model with 78 mixture components.
#>  - Model fitting took 0.04 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>  - Computing 927 x 254 likelihood matrix.
#>  - Likelihood calculations took 0.17 seconds.
#>  - Fitting model with 254 mixture components.
#>  - Model fitting took 0.14 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>  - Computing 927 x 265 likelihood matrix.
#>  - Likelihood calculations took 0.18 seconds.
#>  - Fitting model with 265 mixture components.
#>  - Model fitting took 0.12 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
fit
#> class: PACEFit
#> cell types (6): B_Cell, Endothelial, Fibroblast, Macrophage, T_Cell, Tumour
#> kernels: h_bio = 30 um, h_tech = 5 um | contamination: percell_hc; dispersion: nb1
#> condition: Responder (ResponderPD)
#> pipeline: model -> shrink -> decompose -> drivers
#>   neighbour slopes: 33372 rows (0 at lfsr < 0.05)
```

The interaction term is named from the non-reference level of
`condition_col`:

``` r

fit@params$resp_term
#> [1] "ResponderPD"
```

so `ResponderPD:Tumour` is the extra change in expression per unit of
tumour proximity in the PD arm, on top of the `Tumour` effect shared by
both arms.

## The responder block in the decomposition

[`varianceDecomposition()`](https://ecool50.github.io/PACE/reference/varianceDecomposition.md)
now carries a `Responder spatial %` column beside the `Spatial %` one:

``` r

varianceDecomposition(fit) |>
  group_by(focal) |>
  summarise(`Spatial %`           = round(median(`Spatial %`), 4),
            `Responder spatial %` = round(median(`Responder spatial %`), 4),
            `Spillover %`         = round(median(`Spillover %`), 2),
            .groups = "drop")
#> # A tibble: 6 × 4
#>   focal       `Spatial %` `Responder spatial %` `Spillover %`
#>   <chr>             <dbl>                 <dbl>         <dbl>
#> 1 B_Cell           0.0151                0               6.21
#> 2 Endothelial      0.0021                0               5.3 
#> 3 Fibroblast       0.0007                0.0004          5.78
#> 4 Macrophage       0.0344                0.0028          8.34
#> 5 T_Cell           0.0232                0               6.54
#> 6 Tumour           0.0342                0.0106          1.31
```

The responder block is far smaller than the shared spatial block, which
is the expected shape: most of how a cell responds to its neighbours is
common to both arms, and only a thin slice of it depends on outcome.
Contamination is larger than either, and largest in the sparsely
distributed immune populations.

## Which pairs carry the condition difference

[`pairVariance()`](https://ecool50.github.io/PACE/reference/pairVariance.md)
and
[`plotPairHeatmap()`](https://ecool50.github.io/PACE/reference/plotPairHeatmap.md)
both take `block = "responder"`, which attributes the responder block
across neighbours instead of the shared spatial one:

``` r

pairVariance(fit, block = "responder") |>
  arrange(desc(val)) |>
  head(5)
#> # A tibble: 5 × 3
#>   focal      neighbour       val
#>   <chr>      <chr>         <dbl>
#> 1 Tumour     Fibroblast  0.0278 
#> 2 Macrophage Tumour      0.0164 
#> 3 Tumour     Endothelial 0.0101 
#> 4 Tumour     T_Cell      0.00274
#> 5 Tumour     Macrophage  0.00233

plotPairHeatmap(fit, block = "responder")
```

![](condition_files/figure-html/pairs-responder-1.png)

Macrophages next to tumour cells is among the strongest pairs here,
which is the relationship the accompanying paper reports for the full
cohort. Compare it with the shared spatial block on the same fit:

``` r

plotPairHeatmap(fit, block = "spatial")
```

![](condition_files/figure-html/pairs-spatial-1.png)

## Statistical power, and what this subset cannot show

The condition contrast is between patients, and the outcome-dependent
part of the spatial signal is sparse. Both facts mean it needs the whole
cohort. In this subset, no gene reaches the usual significance
threshold:

``` r

ns <- neighbourSlopes(fit)
c(rows = nrow(ns), calls_lfsr_under_0.05 = sum(ns$lfsr < 0.05, na.rm = TRUE))
#>                  rows calls_lfsr_under_0.05 
#>                 33372                     0
```

[`topDrivers()`](https://ecool50.github.io/PACE/reference/topDrivers.md)
scores every pair and labels each one, and its own verdict for all
thirty here is the right one:

``` r

table(vapply(topDrivers(fit), function(x) x$status, character(1)))
#> 
#> honestly null 
#>            30
```

### What the full cohort gives

Rather than assert the difference, the package ships the macrophage
slopes from the **full** cohort, all 56,274 cells across these same 26
patients, so the two can be put side by side. It is a fitted result, not
something this vignette computes; the `provenance` attribute records
how:

``` r

full <- readRDS(system.file("extdata", "mel_full_cohort_macrophage_slopes.rds",
                            package = "PACE"))
nrow(full)
#> [1] 5562
subset(full, gene == "SPP1" & neighbour == "Tumour" & term == "ResponderPD:Tumour")
#>     gene      focal neighbour               term    estimate   std.error
#> 810 SPP1 Macrophage    Tumour ResponderPD:Tumour -0.06738296 0.006050728
#>     estimate_shrunk   sd_shrunk lfsr
#> 810     -0.06683963 0.006026284    0
```

That is the result the paper reports: in macrophages, the response of
*SPP1* to tumour proximity differs between arms, more negative in
progressive disease. Now the same row from the subset fitted above:

``` r

subset(ns, gene == "SPP1" & focal == "Macrophage" &
           neighbour == "Tumour" & term == "ResponderPD:Tumour")
#>       gene      focal neighbour               term     estimate   std.error
#> 31401 SPP1 Macrophage    Tumour ResponderPD:Tumour -0.004324059 0.006125521
#>       estimate_shrunk sd_shrunk lfsr
#> 31401               0         0    1
```

The estimate has been shrunk to exactly zero. This is worth being
precise about, because it is not a threshold that could be relaxed:
32,445 of the 33,372 responder terms are shrunk to zero here, mash
having concluded that the block carries no signal at this sample size.
No `lfsr` cut recovers *SPP1*, because there is no estimate left to
recover. The unshrunk estimate does keep the right sign, but it is
smaller than its own standard error.

Cropping each section to 15% removes the effect entirely; it only begins
to reappear at around 60% of the cohort, and even there it is less than
half its full size. Crops targeted at the macrophages do recover it, but
only by keeping most of the tissue, because the macrophages are spread
through all of it.

So the pair ranking above is worth reading as a direction, and the
significance column as a warning. For the published result, fit the full
deposit.

The practical lesson generalises beyond this dataset: a condition
analysis is powered by **patients**, not by cells. Adding cells from the
patients you already have sharpens each patient’s estimate, but the
contrast is still a comparison of two small groups.

## Session info

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] dplyr_1.2.1                 ggplot2_4.0.3              
#>  [3] SpatialExperiment_1.22.0    SingleCellExperiment_1.34.0
#>  [5] SummarizedExperiment_1.42.0 Biobase_2.72.0             
#>  [7] GenomicRanges_1.64.0        Seqinfo_1.2.0              
#>  [9] IRanges_2.46.0              S4Vectors_0.50.2           
#> [11] BiocGenerics_0.58.1         generics_0.1.4             
#> [13] MatrixGenerics_1.24.0       matrixStats_1.5.0          
#> [15] PACE_0.99.1                 BiocStyle_2.40.0           
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6        rjson_0.2.23        xfun_0.60          
#>  [4] bslib_0.12.0        lattice_0.22-9      vctrs_0.7.3        
#>  [7] tools_4.6.1         parallel_4.6.1      tibble_3.3.1       
#> [10] pkgconfig_2.0.3     Matrix_1.7-5        SQUAREM_2026.1     
#> [13] RColorBrewer_1.1-3  S7_0.2.2            desc_1.4.3         
#> [16] assertthat_0.2.1    lifecycle_1.0.5     truncnorm_1.0-9    
#> [19] compiler_4.6.1      farver_2.1.2        textshaping_1.0.5  
#> [22] codetools_0.2-20    htmltools_0.5.9     sass_0.4.10        
#> [25] yaml_2.3.12         tidyr_1.3.2         pillar_1.11.1      
#> [28] pkgdown_2.2.1       jquerylib_0.1.4     BiocParallel_1.46.0
#> [31] rmeta_3.0           cachem_1.1.0        DelayedArray_0.38.2
#> [34] dbscan_1.2.6        magick_2.9.1        abind_1.4-8        
#> [37] tidyselect_1.2.1    digest_0.6.39       mvtnorm_1.4-2      
#> [40] purrr_1.2.2         bookdown_0.48       ashr_2.2-63        
#> [43] labeling_0.4.3      fastmap_1.2.0       grid_4.6.1         
#> [46] cli_3.6.6           invgamma_1.2        SparseArray_1.12.2 
#> [49] magrittr_2.0.5      S4Arrays_1.12.0     utf8_1.2.6         
#> [52] withr_3.0.3         scales_1.4.0        rmarkdown_2.32     
#> [55] XVector_0.52.0      otel_0.2.0          ragg_1.5.2         
#> [58] evaluate_1.0.5      knitr_1.52          viridisLite_0.4.3  
#> [61] irlba_2.3.7         rlang_1.3.0         Rcpp_1.1.2         
#> [64] mixsqp_0.3-54       glue_1.8.1          BiocManager_1.30.27
#> [67] jsonlite_2.0.0      plyr_1.8.9          mashr_0.2.79       
#> [70] R6_2.6.1            systemfonts_1.3.2   fs_2.1.0
```
