---
title: "Bisplotti User Guide"
date: "21 February 2023"
package: "bisplotti 0.0.18"
output:
    BiocStyle::html_document:
        highlight: pygments
        toc_float: true
        fig_width: 8
        git_height: 6
vignette: >
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteIndexEntry{Bisplotti User Guide}
    %\VignetteEncoding[utf8]{inputenc}
---



# Bisplotti

`bisplotti` is a package to generate commonly produced plots in DNA methylation
sequencing analyses. It creates plots from standard Bioconductor formats (i.e.,
`GRanges`) and the commonly used `BSseq` format.

# Quick Start

## Installing

A development version is available on GitHub and can be installed via:


```r
if (!requireNamespace("BiocManager", quietly=TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("huishenlab/bisplotti")
```


```r
library(bisplotti)
```

# Create Plots

## Epiread

To create the lollipop epiread plots, two commands are needed. First, use the
GRanges output from `biscuiteer::readEpibed()` as input to `tabulateEpibed`,
which turns the GRanges into a convenient matrix format. The second command,
`plotEpiread` uses the matrix output from `tabulateEpibed` as the input to
create the lollipop plot. `epistateCaller()` can take the output of
`tabulateEpibed()` to cluster the epireads on both CpG and GpC methylation.
Note, `epistateCaller()` is under development and is not currently suggested
for use in publication level analyses.


```r
epibed.nome     <- system.file("extdata", "hct116.nome.epibed.gz", package="biscuiteer")
epibed.nome.gr  <- readEpibed(epibed = epibed.nome, genome = "hg19", chr = "chr1")
epibed.tab.nome <- tabulateEpibed(epibed.nome.gr)
plotEpiread(epibed.tab.nome$gc_table)
```

```
## Error in plotEpiread(epibed.tab.nome$gc_table): could not find function "plotEpiread"
```

```r
#epistateCaller(epibed.tab.nome)
```

## Multiscale

The multiscale plot is based on Figures [2](https://www.nature.com/articles/s41588-018-0073-4/figures/2) 
and [6](https://www.nature.com/articles/s41588-018-0073-4/figures/6) of 
[Zhou et al., Nature Genetics, 2018](https://www.nature.com/articles/s41588-018-0073-4).
It shows the average methylation levels across varying size genomic windows.
The example shows a small portion of chromosome 16 for window sizes running
from 1 Mb to 10 Mb.


```r
files.loc <- system.file("extdata", package="bisplotti")

files <- lapply(
    list.files(files.loc, pattern="Heyn_2012_100yr", full.names=TRUE), function(x) {
        return(rtracklayer::import(x, format="bedGraph"))
    }
)
cnames <- list.files(files.loc, pattern="Heyn_2012_100yr")
cnames <- gsub("Heyn_2012_100yr", "100yr", cnames)
cnames <- gsub(".bed.gz", "", cnames)
names(files) <- cnames

files.grl <- as(files, "GRangesList")

multiscaleMethylationPlot(files.grl)
```

![plot of chunk multiscale](figure/multiscale-1.png)

## 1D Methylation Level Density

The 1D methylation level density plot shows the density of the methylation
levels for all CpGs in your dataset. An example of creating this plot is:


```r
bed <- system.file("extdata", "MCF7_Cunha_chr11p15.bed.gz", package="biscuiteer")
vcf <- system.file("extdata", "MCF7_Cunha_header_only.vcf.gz", package="biscuiteer")
bisc <- readBiscuit(BEDfile = bed, VCFfile = vcf, merged = FALSE)

meth1DDensity(bisc)
```

![plot of chunk meth1d](figure/meth1d-1.png)
A matrix of methylation levels (beta values) can also be used as input to
`meth1DDensity()`. This method assumes the column names of your matrix are the
samples in the matrix. The row names can either be NULL or the CpG loci/probes.

## 2D Methylation Level Density

The 2D methylation level density plot compares the density of average
methylation levels in provided bins across two samples. An example of creating
this plot is:


```r
orig_bed <- system.file("extdata", "MCF7_Cunha_chr11p15.bed.gz", package="biscuiteer")
orig_vcf <- system.file("extdata", "MCF7_Cunha_header_only.vcf.gz", package="biscuiteer")
bisc1    <- readBiscuit(BEDfile = orig_bed, VCFfile = orig_vcf, merged = FALSE)

shuf_bed <- system.file("extdata", "MCF7_Cunha_chr11p15_shuffled.bed.gz", package="biscuiteer")
shuf_vcf <- system.file("extdata", "MCF7_Cunha_shuffled_header_only.vcf.gz", package="biscuiteer")
bisc2    <- readBiscuit(BEDfile = shuf_bed, VCFfile = shuf_vcf, merged = FALSE)

comb <- unionize(bisc1, bisc2)

meth2DDensity(comb, chr="chr11", sample_1 = "MCF7_Cunha", sample_2 = "MCF7_Cunha_shuffled")
```

![plot of chunk meth2d](figure/meth2d-1.png)

A matrix of methylation levels (beta values) can also be used as input to
`meth2DDensity()`. This method assumes the column names of your matrix are the
samples in the matrix. The row names must be the CpG loci/probes (i.e., the
result of `rownames(mat) <- as.character(granges(gr))`. In either input case,
there must be at least two samples in your input object.

# Session Info


```r
sessionInfo()
```

```
## R version 4.2.1 (2022-06-23)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur ... 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] bisplotti_0.0.18            RColorBrewer_1.1-3         
##  [3] ggplot2_3.4.1               biscuiteer_1.13.0          
##  [5] bsseq_1.34.0                SummarizedExperiment_1.28.0
##  [7] Biobase_2.58.0              MatrixGenerics_1.10.0      
##  [9] matrixStats_0.63.0          GenomicRanges_1.50.2       
## [11] GenomeInfoDb_1.34.9         IRanges_2.32.0             
## [13] S4Vectors_0.36.1            biscuiteerData_1.12.0      
## [15] ExperimentHub_2.6.0         AnnotationHub_3.6.0        
## [17] BiocFileCache_2.6.1         dbplyr_2.3.0               
## [19] BiocGenerics_0.44.0        
## 
## loaded via a namespace (and not attached):
##   [1] plyr_1.8.8                               
##   [2] splines_4.2.1                            
##   [3] BiocParallel_1.32.5                      
##   [4] listenv_0.9.0                            
##   [5] CGHcall_2.60.0                           
##   [6] digest_0.6.31                            
##   [7] foreach_1.5.2                            
##   [8] htmltools_0.5.4                          
##   [9] viridis_0.6.2                            
##  [10] GO.db_3.16.0                             
##  [11] fansi_1.0.4                              
##  [12] magrittr_2.0.3                           
##  [13] memoise_2.0.1                            
##  [14] BSgenome_1.66.3                          
##  [15] tzdb_0.3.0                               
##  [16] limma_3.54.1                             
##  [17] readr_2.1.4                              
##  [18] globals_0.16.2                           
##  [19] Biostrings_2.66.0                        
##  [20] R.utils_2.12.2                           
##  [21] prettyunits_1.1.1                        
##  [22] colorspace_2.1-0                         
##  [23] blob_1.2.3                               
##  [24] rappdirs_0.3.3                           
##  [25] xfun_0.37                                
##  [26] dplyr_1.1.0                              
##  [27] crayon_1.5.2                             
##  [28] RCurl_1.98-1.10                          
##  [29] org.Mm.eg.db_3.16.0                      
##  [30] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2  
##  [31] graph_1.76.0                             
##  [32] annotatr_1.24.0                          
##  [33] Mus.musculus_1.3.1                       
##  [34] impute_1.72.3                            
##  [35] iterators_1.0.14                         
##  [36] VariantAnnotation_1.44.1                 
##  [37] glue_1.6.2                               
##  [38] gtable_0.3.1                             
##  [39] zlibbioc_1.44.0                          
##  [40] XVector_0.38.0                           
##  [41] DelayedArray_0.24.0                      
##  [42] dmrseq_1.18.0                            
##  [43] Rhdf5lib_1.20.0                          
##  [44] future.apply_1.10.0                      
##  [45] HDF5Array_1.26.0                         
##  [46] scales_1.2.1                             
##  [47] rngtools_1.5.2                           
##  [48] DBI_1.1.3                                
##  [49] CGHbase_1.58.0                           
##  [50] Rcpp_1.0.10                              
##  [51] viridisLite_0.4.1                        
##  [52] xtable_1.8-4                             
##  [53] progress_1.2.2                           
##  [54] bumphunter_1.40.0                        
##  [55] bit_4.0.5                                
##  [56] OrganismDbi_1.40.0                       
##  [57] httr_1.4.4                               
##  [58] ellipsis_0.3.2                           
##  [59] farver_2.1.1                             
##  [60] pkgconfig_2.0.3                          
##  [61] XML_3.99-0.13                            
##  [62] R.methodsS3_1.8.2                        
##  [63] locfit_1.5-9.7                           
##  [64] utf8_1.2.3                               
##  [65] DNAcopy_1.72.3                           
##  [66] labeling_0.4.2                           
##  [67] reshape2_1.4.4                           
##  [68] tidyselect_1.2.0                         
##  [69] rlang_1.0.6                              
##  [70] later_1.3.0                              
##  [71] AnnotationDbi_1.60.0                     
##  [72] munsell_0.5.0                            
##  [73] BiocVersion_3.16.0                       
##  [74] tools_4.2.1                              
##  [75] cachem_1.0.6                             
##  [76] cli_3.6.0                                
##  [77] generics_0.1.3                           
##  [78] RSQLite_2.3.0                            
##  [79] evaluate_0.20                            
##  [80] stringr_1.5.0                            
##  [81] fastmap_1.1.0                            
##  [82] yaml_2.3.7                               
##  [83] outliers_0.15                            
##  [84] org.Hs.eg.db_3.16.0                      
##  [85] knitr_1.42                               
##  [86] bit64_4.0.5                              
##  [87] purrr_1.0.1                              
##  [88] KEGGREST_1.38.0                          
##  [89] doRNG_1.8.6                              
##  [90] nlme_3.1-162                             
##  [91] RBGL_1.74.0                              
##  [92] future_1.31.0                            
##  [93] sparseMatrixStats_1.10.0                 
##  [94] mime_0.12                                
##  [95] R.oo_1.25.0                              
##  [96] xml2_1.3.3                               
##  [97] biomaRt_2.54.0                           
##  [98] BiocStyle_2.26.0                         
##  [99] compiler_4.2.1                           
## [100] filelock_1.0.2                           
## [101] curl_5.0.0                               
## [102] png_0.1-8                                
## [103] interactiveDisplayBase_1.36.0            
## [104] marray_1.76.0                            
## [105] tibble_3.1.8                             
## [106] Homo.sapiens_1.3.1                       
## [107] stringi_1.7.12                           
## [108] QDNAseq_1.34.0                           
## [109] highr_0.10                               
## [110] GenomicFeatures_1.50.4                   
## [111] lattice_0.20-45                          
## [112] Matrix_1.5-3                             
## [113] permute_0.9-7                            
## [114] vctrs_0.5.2                              
## [115] pillar_1.8.1                             
## [116] lifecycle_1.0.3                          
## [117] rhdf5filters_1.10.0                      
## [118] BiocManager_1.30.19                      
## [119] cowplot_1.1.1                            
## [120] data.table_1.14.8                        
## [121] bitops_1.0-7                             
## [122] TxDb.Mmusculus.UCSC.mm10.knownGene_3.10.0
## [123] httpuv_1.6.9                             
## [124] rtracklayer_1.58.0                       
## [125] R6_2.5.1                                 
## [126] BiocIO_1.8.0                             
## [127] promises_1.2.0.1                         
## [128] gridExtra_2.3                            
## [129] KernSmooth_2.23-20                       
## [130] parallelly_1.34.0                        
## [131] codetools_0.2-19                         
## [132] MASS_7.3-58.2                            
## [133] gtools_3.9.4                             
## [134] assertthat_0.2.1                         
## [135] qualV_0.3-4                              
## [136] rhdf5_2.42.0                             
## [137] rjson_0.2.21                             
## [138] withr_2.5.0                              
## [139] regioneR_1.30.0                          
## [140] GenomicAlignments_1.34.0                 
## [141] Rsamtools_2.14.0                         
## [142] GenomeInfoDbData_1.2.9                   
## [143] parallel_4.2.1                           
## [144] hms_1.1.2                                
## [145] grid_4.2.1                               
## [146] rmarkdown_2.20                           
## [147] DelayedMatrixStats_1.20.0                
## [148] shiny_1.7.4                              
## [149] restfulr_0.0.15
```

