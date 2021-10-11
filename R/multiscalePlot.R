#' Create the multiscale plot from a list of GRanges
#' 
#' @param grl GRangesList of a single set of multiscale BED files
#' @param colors RColorBrewer color scheme to use (default: YlGnBu)
#' @param na_color Color to give NA values (default: darkgray)
#' @param what Column to plot (default: score)
#' @param chr.to.plot Chromosome to plot, all if NULL (default: NULL)
#' 
#' @import GenomicRanges
#' @import ggplot2
#' @import viridis
#' 
#' @export
#' 
#' @example
#' 
multiscaleMethylationPlot <- function(grl,
                                      colors="YlGnBu",
                                      na_color="darkgray",
                                      what="score",
                                      chr.to.plot = NULL) {

    # Check inputs are in correct form
    if (!is(grl, "GRangesList")) stop("Input needs to be a GRangesList object.")
    stopifnot(what %in% names(mcols(grl[[1]])))

    # Pull out chromosome to plot, if requested
    if (!is.null(chr.to.plot)) {
        grl <- lapply(grl, function(x) { return (keepSeqlevels(x, chr.to.plot, "coarse"))})
        grl <- as(grl, "GRangesList")
    }

    # Turn GRangesList into matrix-like object for easy plotting
    score.matrix <- matrix(
        0,
        nrow = sort(unlist(lapply(grl, length)), decreasing = TRUE)[1],
        ncol = length(grl)
    )

    long.mat <- order(unlist(lapply(grl, length)), decreasing = TRUE)[1]
    long.mat <- grl[[long.mat]]
    rownames(score.matrix) <- as.character(granges(long.mat))
    colnames(score.matrix) <- names(grl)

    score.matrix <- as.data.frame(
        do.call(
            cbind,
            lapply(
                1:length(grl),
                function(x) {
                    meta <- mcols(grl[[x]])
                    rownames(meta) <- as.character((granges(grl[[x]])))
                    score.matrix.sub <- as.matrix(score.matrix[,names(grl)[x]])
                    score.matrix.sub[rownames(meta),] <- meta$score
                    return(score.matrix.sub)
                }
            )
        )
    )
    colnames(score.matrix) <- names(grl)

    score.matrix <- as.data.frame(t(score.matrix))
    score.matrix$res <- rownames(score.matrix)
    score.matrix.melt <- reshape2::melt(score.matrix, id.vars = "res")
    score.matrix.melt$res <- levels(factor(score.matrix.melt$res))

    # Create plot
    g <- ggplot(score.matrix.melt, aes(x = variable, y = res, fill = value)) +
        geom_raster(interpolate = TRUE) +
        #scale_fill_gradientn(colours = colorRampPalette(c("blue","yellow"))(20), na.) +
        scale_fill_distiller(palette = colors, direction = -1, na.value = na_color,
                             name = "Methylation level") +
        scale_y_discrete(name = "Resolution",
                           breaks = c("New_chr16_4.0", "New_chr16_5.0",
                                      "New_chr16_6.0", "New_chr16_7.0"),
                           labels = c("10kb", "100kb", "1Mb", "10Mb")) +
        theme_bw(12) +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            #axis.title.y = element_blank(),
            #axis.text.y = element_blank(),
            #axis.ticks.y = element_blank(),
            #axis.text = element_blank(),
            #axis.ticks = element_blank(),
            #legend.title = element_blank(),
        )
  
  return(g)
}

# Get the mean multiscale data
# setwd("~/Downloads")
# means <- lapply(
#     list.files("multiscale_analysis/stats", pattern="cross_sample_mean", full.names=TRUE), function(x) {
#         return(rtracklayer::import(x, format = "bedGraph"))
#     }
# )
# cnames <- list.files("multiscale_analysis/stats", pattern="cross_sample_mean")
# cnames <- gsub("cross_sample_mean", "mean", cnames)
# cnames <- gsub(".bed.gz", "", cnames)
# names(means) <- cnames
# 
# # Only want the first 35 Mb for comparison with Nature Genetics Paper
# chr16_0_35Mb <- GRanges(Rle(c("chr16"), c(1)), IRanges(1, width = 35000000))
# 
# # Single samples
# sampl_100 <- lapply(
#     list.files("multiscale_analysis/means", pattern="Heyn_2012_Human_CD4T_100yo", full.names=TRUE), function(x) {
#         return(rtracklayer::import(x, format = "bedGraph", which = chr16_0_35Mb))
#     }
# )
# cnames <- list.files("multiscale_analysis/means", pattern="Heyn_2012_Human_CD4T_100yo")
# cnames <- gsub("Heyn_2012_Human_CD4T_100yo", "100yr", cnames)
# cnames <- gsub(".bed.gz", "", cnames)
# names(sampl_100) <- cnames
# 
# sampl_New <- lapply(
#     list.files("multiscale_analysis/means", pattern="Heyn_2012_Human_CD4T_Newborn", full.names=TRUE), function(x) {
#         return(rtracklayer::import(x, format = "bedGraph", which = chr16_0_35Mb))
#     }
# )
# cnames <- list.files("multiscale_analysis/means", pattern="Heyn_2012_Human_CD4T_Newborn")
# cnames <- gsub("Heyn_2012_Human_CD4T_Newborn", "New", cnames)
# cnames <- gsub(".bed.gz", "", cnames)
# names(sampl_New) <- cnames
# 
# # Turn list of GRanges into GRangesList
# test.data.means <- as(means, "GRangesList")
# test.data.sampl_100 <- as(sampl_100, "GRangesList")
# test.data.sampl_New <- as(sampl_New, "GRangesList")
# 
# # Create multiscale methylation plots
# pdf(file="test/test_means.pdf", width = 8, height = 3)
# multiscaleMethylationPlot(test.data.means, chr.to.plot = "chr16", colors = "YlGnBu")
# dev.off()
# 
# pdf(file="test/test_sampl_100.pdf", width = 8, height = 3)
# multiscaleMethylationPlot(test.data.sampl_100, chr.to.plot = NULL, colors = "YlGnBu")
# dev.off()
# 
# pdf(file="test/test_sampl_New.pdf", width = 8, height = 3)
# multiscaleMethylationPlot(test.data.sampl_New, chr.to.plot = NULL, colors = "YlGnBu")
# dev.off()
