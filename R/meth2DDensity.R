#' Create a 2D density plot of methylation levels for two samples
#'
#' @param betas BSseq or matrix of methylation values, must contain at least two samples
#' @param sample_1 Sample name to plot on x-axis, sample at index 1 if NULL (default: NULL)
#' @param sample_2 Sample name to plot on y-axis, sample at index 2 if NULL (default: NULL)
#' @param which Integer or GRanges. If Integer, find bins using tileGenome() with Integer tilewidth, otherwise use
#' GRanges. If GRanges, ignores genome and chr (default: 10000)
#' @param genome Genome data is from, currently only works for hg19, hg38, mm9, and mm10 (default: "hg38")
#' @param chr Pull out a single chromosome for bins, use all chromosomes if NULL (default: NULL)
#' @param palette RColorBrewer palette to create plot with (default: "PuBu")
#' @param make_log Make z-axis for density be log10 (default: FALSE)
#'
#' @return a ggplot2 object of the density plot
#'
#' @import GenomicRanges
#' @import biscuiteer
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom cowplot theme_cowplot
#'
#' @export
#'
#' @examples
#'
#' orig_bed <- system.file("extdata", "MCF7_Cunha_chr11p15.bed.gz", package="biscuiteer")
#' orig_vcf <- system.file("extdata", "MCF7_Cunha_header_only.vcf.gz", package="biscuiteer")
#' bisc1    <- readBiscuit(BEDfile = orig_bed, VCFfile = orig_vcf, merged = FALSE)
#' 
#' shuf_bed <- system.file("extdata", "MCF7_Cunha_chr11p15_shuffled.bed.gz", package="biscuiteer")
#' shuf_vcf <- system.file("extdata", "MCF7_Cunha_shuffled_header_only.vcf.gz", package="biscuiteer")
#' bisc2    <- readBiscuit(BEDfile = shuf_bed, VCFfile = shuf_vcf, merged = FALSE)
#'
#' comb <- unionize(bisc1, bisc2)
#'
#' meth2DDensity(comb, chr="chr11", sample_1 = "MCF7_Cunha", sample_2 = "MCF7_Cunha_shuffled")
#'
meth2DDensity <- function(betas,
                          sample_1 = NULL,
                          sample_2 = NULL,
                          which = 10000,
                          genome = "hg38",
                          chr = NULL,
                          palette = "PuBu",
                          make_log = FALSE) {

    # Check inputs are in correct form
    if (!is(betas, "BSseq") & !is(betas, "matrix")) stop("betas needs to be a BSseq or a matrix object.")
    if (!ncol(betas) >= 2) stop("betas must have at least two samples")
    if (!is(which, "numeric") & !is(which, "GRanges")) stop("which must be an integer or a GRanges object.")
    if (!(genome %in% c("hg19", "hg38", "mm10")) & !is(which, "GRanges")) {
        stop("genome unrecognized. Please provide a GRanges with your bins.")
    }
    if (!is.logical(make_log)) stop("make_log must be TRUE or FALSE")

    if (is(which, "numeric")) {
        gen.sl <- getSeqLengths(genome=genome, chr=chr)
        bins <- tileGenome(seqlengths = gen.sl, tilewidth = which, cut.last.tile.in.chrom = TRUE)
        seqlengths(bins) <- gen.sl
    } else {
        bins <- which
    }

    ## Setup betas matrix
    if (is(betas, "BSseq")) {
        meth <- getMeth(betas, regions = bins, type="raw", what="perRegion")
    } else {
        if (!is(rownames(betas), "character")) stop("betas row names must be loci positions (chrA:1234)")
        tmp <- as(betas, "GRanges")
        meth <- averagePerBinBSseqLike(tmp, bins)
    }

    if (!is.null(sample_1)) {
        if (!(sample_1 %in% colnames(meth))) {
            idx1 <- 1
            message("Cannot find ", sample_1, ". Setting x-axis sample to ", colnames(meth)[idx1])
        } else {
            idx1 <- match(sample_1, colnames(meth))
        }
    } else {
        idx1 <- 1
    }

    if (!is.null(sample_2)) {
        if (!(sample_2 %in% colnames(meth))) {
            idx2 <- 2
            message("Cannot find ", sample_2, ". Setting y-axis sample to ", colnames(meth)[idx2])
        } else {
            idx2 <- match(sample_2, colnames(meth))
        }
    } else {
        idx2 <- 2
    }

    df <- as.data.frame(meth)
    df <- df[complete.cases(df),] # drop rows that don't have data in both samples
    s1 <- colnames(df)[idx1]
    s2 <- colnames(df)[idx2]

    g <- ggplot(mapping=aes(x=df[[s1]], y=df[[s2]]) )

    if (make_log) {
        # Find minimum color to set as NA value for log scale
        cols <- brewer.pal(n = 8, name = palette)

        g <- g +
            stat_density_2d(aes(fill = log10(..density..)), geom = "raster", contour = FALSE, na.rm=TRUE) +
            scale_fill_distiller(palette = palette, direction = -1, limits=c(-3, NA), na.value = cols[length(cols)])
    } else {
        g <- g +
            stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, na.rm=TRUE) +
            scale_fill_distiller(palette = palette, direction = -1)
    }

    g <- g +
        scale_x_continuous(breaks=seq(0, 1, 0.2), limits = c(0, 1)) +
        scale_y_continuous(breaks=seq(0, 1, 0.2), limits = c(0, 1)) +
        xlab(s1) +
        ylab(s2) +
        cowplot::theme_cowplot()

    return(g)
}
