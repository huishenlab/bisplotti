#' Merge metadata columns of two GRanges
#'
#' Merges mcols(gr2) into mcols(gr1). NA's are set for rows that occur in one GRanges but not the other. Once merged,
#' the output GRanges will be sorted.
#'
#' @param gr1 GRanges to merge metadata into
#' @param gr2 GRanges with data to be merged
#' @param append.gr1 String to append to gr1 columns that share names with gr2 (default: "_1")
#' @param append.gr2 String to append to gr2 columns that share names with gr1 (default: "_2")
#'
#' @return GRanges
#'
#' @import GenomicRanges
#'
#' @export
#'
#' @examples
#'
#' gr1 <- GRanges(seqnames = "chr2", ranges = IRanges(103, 106), strand = "+", score = 5L, GC = 0.45)
#' gr2 <- GRanges(
#'     seqnames = c("chr1", "chr2"),
#'     ranges = IRanges(c(107, 113), width = 3),
#'     strand = c("+", "-"),
#'     score = 3:4, GC = c(0.3, 0.5))
#' gr3 <- mergeMCols(gr1, gr2)
#'
mergeMCols <- function(gr1, gr2, append.gr1 = "_1", append.gr2 = "_2") {
    # Check inputs
    if (!is(gr1, "GenomicRanges")) stop("gr1 must be a GenomicRanges object")
    if (!is(gr2, "GenomicRanges")) stop("gr2 must be a GenomicRanges object")
    if (!is(append.gr1, "character") & !length(append.gr1) == 1) stop("append.gr1 must be a character object of length 1")
    if (!is(append.gr2, "character") & !length(append.gr2) == 1) stop("append.gr2 must be a character object of length 1")

    # Check if any column names are the same, tack append.gr1 and append.gr2 to any shared column names
    # Note, this isn't very "smart" in that it won't gracefully handle names that are already shared
    #     (i.e., name_1 in both gr1 and gr2 will become name_1_1 and name_1_2 after combining)
    gr1.names <- colnames(mcols(gr1))
    gr2.names <- colnames(mcols(gr2))

    gr1.names.new <- unlist(lapply(gr1.names, function(x) { ifelse(x %in% gr2.names, paste0(x, append.gr1), x) }))
    gr2.names.new <- unlist(lapply(gr2.names, function(x) { ifelse(x %in% gr1.names, paste0(x, append.gr2), x) }))
    colnames(mcols(gr1)) <- gr1.names.new
    colnames(mcols(gr2)) <- gr2.names.new

    # If the granges of the two inputs are the same, a simple cbind should suffice, however this also requires
    #     the seqinfo is the same for both as well.
    # It's overall easier to convert to data.frames, merge, then convert back
    df1 <- as.data.frame(gr1)
    df2 <- as.data.frame(gr2)
    df3 <- merge(df1, df2,
                 by = c("seqnames", "start", "end", "strand", "width"),
                 all = TRUE, suffixes = c(append.gr1, append.gr2))
    gr1 <- makeGRangesFromDataFrame(df3, keep.extra.columns = TRUE)

    return(sort(gr1))
}

#' Get the seqlengths of a chromosome or all chromosomes
#'
#' @param genome Genome to retrieve, currently only hg38, hg19, mm9, mm10 (default: hg38)
#' @param chr Chromosome to retrieve length of, all chromosomes if NULL (default: NULL)
#'
#' @return seqlengths of genome
#'
#' @import GenomicRanges
#'
#' @export
#'
#' @examples
#'
#' gen.sl <- getSeqLengths("hg38", "chr11")
#'
getSeqLengths <- function(genome = c("hg38", "hg19", "mm9", "mm10"), chr = NULL) {
    # Eventually we should support arbitrary genomes
    genome <- match.arg(genome)

    # Check if the genome used exists in what is currently supported, stopping if not
    if (!genome %in% c("hg38", "hg19", "mm9", "mm10")) stop("Only human and mouse are supported for the time being.")

    # Import
    genome.gr <- switch(genome,
                        hg19 = data("hg19.gr", package = "bisplotti"),
                        hg38 = data("hg38.gr", package = "bisplotti"),
                        mm9  = data("mm9.gr" , package = "bisplotti"),
                        mm10 = data("mm10.gr", package = "bisplotti"))

    # Make sure that the chromosome specified exists in the seqlevels
    if (!is.null(chr)) {
        if (!chr %in% seqlevels(get(genome.gr))) stop("Desired chromosome is not found in the seqlevels of ", genome)
    }

    # Get the seqlengths
    if (is.null(chr)) {
        sl <- seqlengths(get(genome.gr))
    } else {
        sl <- seqlengths(get(genome.gr))[chr]
    }

    return(sl)
}

#' Find the average values of x in provided bins
#'
#' Uses a method similar to the GenomicRanges tutorial for finding binned averages. Note, this method may not return the
#' same results as averagePerBinBSseqLike().
#'
#' @param x GRanges with values to find average of
#' @param bins GRanges with bins to average over
#'
#' @return matrix
#'
#' @import GenomicRanges
#'
#' @export
#'
#' @examples
#'
#' gr1 <- GRanges(seqnames = "chr2", ranges = IRanges(103, 106), strand = "+", score = 5L, GC = 0.45)
#' gr2 <- GRanges(
#'     seqnames = c("chr1", "chr2"),
#'     ranges = IRanges(c(107, 113), width = 3),
#'     strand = c("+", "-"),
#'     score = 3:4, GC = c(0.3, 0.5))
#' gr3 <- mergeMCols(gr1, gr2)
#'
#' bins <- GRanges(seqnames=c("chr1","chr2"), ranges=IRanges(c(100, 100), width=15))
#'
#' averagePerBin(gr3, bins)
#'
averagePerBin <- function(x, bins) {
    # Check inputs
    if (!is(bins, "GenomicRanges")) stop("'bins' must be a GenomicRanges object")
    if (!is(x, "GenomicRanges")) {
        if (is(x, "BSseq")) {
            message("Input a BSseq object, returning bsseq::getMeth(x, regions=bins, what=\"perRegion\", type=\"raw\"")
            return(getMeth(x, regions=bins, what="perRegion", type="raw"))
        } else {
            stop("'x' must be a GenomicRanges object")
        }
    }

    tmp <- lapply(colnames(mcols(x)),
                  FUN = function(m) {
                      t.score <- GenomicRanges::coverage(x, weight=m)
                      t.avg <- GenomicRanges::binnedAverage(bins, t.score, m)
                      return(t.avg)
                  })

    t <- Reduce(mergeMCols, tmp)

    return(t)
}

#' Find the average values of x in provided bins
#'
#' Uses a method similar to BSseq::getMeth(bsseq, regions=bins, what="perRegion"). Note, this method may not return the
#' same results as averagePerBin().
#'
#' @param x GRanges with values to find average of
#' @param bins GRanges with bins to average over
#'
#' @return matrix
#'
#' @import GenomicRanges
#' @importFrom DelayedMatrixStats colMeans2
#'
#' @export
#'
#' @examples
#'
#' gr1 <- GRanges(seqnames = "chr2", ranges = IRanges(103, 106), strand = "+", score = 5L, GC = 0.45)
#' gr2 <- GRanges(
#'     seqnames = c("chr1", "chr2"),
#'     ranges = IRanges(c(107, 113), width = 3),
#'     strand = c("+", "-"),
#'     score = 3:4, GC = c(0.3, 0.5))
#' gr3 <- mergeMCols(gr1, gr2)
#'
#' bins <- GRanges(seqnames=c("chr1","chr2", "chr3"), ranges=IRanges(c(100, 100, 100), width=15))
#'
#' averagePerBinBSseqLike(gr3, bins)
#'
averagePerBinBSseqLike <- function(x, bins) {
    # Check inputs
    if (!is(bins, "GenomicRanges")) stop("'bins' must be a GenomicRanges object")
    if (!is(x, "GenomicRanges")) {
        if (is(x, "BSseq")) {
            message("Input a BSseq object, returning bsseq::getMeth(x, regions=bins, what=\"perRegion\", type=\"raw\"")
            return(getMeth(x, regions=bins, what="perRegion", type="raw"))
        } else {
            stop("'x' must be a GenomicRanges object")
        }
    }

    ol   <- findOverlaps(x, bins)
    meth <- as.matrix(mcols(x)[queryHits(ol), , drop = FALSE])

    out <- lapply(split(meth, subjectHits(ol)), matrix, ncol=ncol(meth))
    out <- do.call(rbind, lapply(out, colMeans2, na.rm=TRUE))

    outMatrix <- matrix(NA, ncol=ncol(mcols(x)), nrow = length(bins))
    outMatrix[as.integer(rownames(out)), ] <- out
    colnames(outMatrix) <- colnames(mcols(x))

    return(outMatrix)
}
