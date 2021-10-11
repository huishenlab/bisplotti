#' Create the multiscale plot from a list of GRanges
#' 
#' @param grl GRangesList of a single set of multiscale BED files
#' @param colors RColorBrewer color scheme to use (default: YlGnBu)
#' @param na_color Color to give NA values (default: darkgray)
#' @param what Column to plot (default: score)
#' @param chr.to.plot Chromosome to plot, all if NULL (default: NULL)
#' 
#' @return a ggplot2 object of the multiscale plot
#' 
#' @import GenomicRanges
#' @import ggplot2
#' @import viridis
#' 
#' @export
#' 
#' @examples
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

