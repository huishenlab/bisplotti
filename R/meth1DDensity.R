#' Create a 1D density plot of methylation level
#'
#' @param betas BSseq or matrix of methylation levels (assumes matrix colnames are the sample names)
#' @param line.color Color of plotted line (default: "red")
#' @param ylabel Y-axis label (default: "Sample")
#'
#' @return a ggplot2 object of the density plot
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom cowplot theme_cowplot
#'
#' @export
#'
#' @examples
#'
#' orig_bed <- system.file("extdata", "MCF7_Cunha_chr11p15.bed.gz",
#'                         package="biscuiteer")
#' orig_vcf <- system.file("extdata", "MCF7_Cunha_header_only.vcf.gz",
#'                         package="biscuiteer")
#' bisc <- readBiscuit(BEDfile = orig_bed, VCFfile = orig_vcf,
#'                     merged = FALSE)
#' meth1DDensity(bisc)
#'
meth1DDensity <- function(betas,
                          line.color="red",
                          ylabel="Sample") {

    # Check inputs are in correct form
    if (!is(betas, "BSseq") & !is(betas, "matrix")) stop("Input needs to be a BSseq or a matrix object.")

    if (is(betas, "BSseq")) {
        meth <- getMeth(betas, type="raw")
    } else {
        meth <- betas
    }
    betas.t <- as.data.frame(t(meth))
    betas.t$sample <- rownames(betas.t)

    betas.melted <- reshape2::melt(betas.t, id.vars = "sample")

    g <- ggplot(betas.melted, aes(x = value, group = sample)) +
        geom_density(fill = NA, color = line.color) +
        xlab("Methylation Level") +
        scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(-0.1,1.1)) +
        ylab(ylabel) +
        cowplot::theme_cowplot() +
        theme(
            legend.position = "none"
        )

    return(g)
}
