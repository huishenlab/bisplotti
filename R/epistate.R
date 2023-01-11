#' Tabulate to a read or fragment level matrix for CGH or GCH context methylation.
#'
#' The input for this function is the GRanges output from biscuiteer::readEpibed().
#' It will auto-detect whether this a BS-seq or a NOMe-seq run and tabulate accordingly.
#'
#' @param gr The epibed GRanges object from readEpibed()
#' @param region Either a GRanges of regions to subset to or explicit region (e.g., chr6:1555-1900)
#' @param include_empty_reads Whether to include reads that contain no methylated (or SNP) sites (default: FALSE)
#' @param include_snps Whether to include SNPs or not, only matters for epiBEDs from BISCUIT v1.0 or 1.1 (default: TRUE)
#'
#' @return A matrix or list of matrices
#'
#' @export
#'
#' @import GenomicRanges
#' @import stringr
#'
#' @examples
#'
#' epibed.nome <- system.file("extdata", "hct116.nome.epiread.gz",
#'                            package="biscuiteer")
#' epibed.nome.gr <- readEpibed(epibed = epibed.nome,
#'                              genome = "hg19", chr = "chr1")
#' epibed.tab.nome <- tabulateEpibed(epibed.nome.gr)
#'
tabulateEpibed <- function(gr,
                           region = NULL,
                           include_empty_reads = FALSE,
                           include_snps = TRUE) {

    # check if a GRanges
    stopifnot(is(gr, "GRanges"))

    # tabulate strings
    #const char skip_epi = '-';
    #const char skip_ins = 'i';
    #const char skip_del = 'd';
    #const char filtered = 'F';
    #const char ignored  = 'x';
    #const char deletion = 'D';
    #const char softclip = 'P';
    #const char methylat = 'M';
    #const char unmethyl = 'U';
    #const char open_acc = 'O';
    #const char shut_acc = 'S';
    #const char ambig_ga = 'R';
    #const char ambig_ct = 'Y';

    # we need to build up a table of CGH and GCH
    if (!all(is.na(gr$CG_decode))) {
        cg_table <- .tabulateRLE(gr, col = "CG_decode", include_snps = include_snps, type = "cg")
        cg_table <- .filterToRegion(cg_table, region = region)
    } else {
        cg_table <- matrix(data = NA, nrow = length(gr$readname), ncol = 1)
        rownames(cg_table) <- gr$readname
        colnames(cg_table) <- "cg_missing"
    }

    if (!all(is.na(gr$GC_decode))) {
        gc_table <- .tabulateRLE(gr, col = "GC_decode", include_snps = include_snps, type = "gc")
        gc_table <- .filterToRegion(gc_table, region = region)
    } else {
        gc_table <- matrix(data = NA, nrow = length(gr$readname), ncol = 1)
        rownames(gc_table) <- gr$readname
        colnames(gc_table) <- "gc_missing"
    }

    if (!all(is.na(gr$VAR_decode))) {
        vr_table <- .tabulateRLE(gr, col = "VAR_decode", include_snps = include_snps, type = "variant")
        vr_table <- .filterToRegion(vr_table, region = region)
    } else {
        vr_table <- matrix(data = NA, nrow = length(gr$readname), ncol = 1)
        rownames(vr_table) <- gr$readname
        colnames(vr_table) <- "variant_missing"
    }

    # Match row order up to the CG table
    reorder_idx_gc <- match(rownames(gc_table), rownames(cg_table))
    reorder_idx_vr <- match(rownames(vr_table), rownames(cg_table))
    gc_table <- gc_table[reorder_idx_gc, , drop = FALSE]
    vr_table <- vr_table[reorder_idx_vr, , drop = FALSE]

    # Remove any rows that can't be found in all tables
    if (!include_empty_reads) {
        drop_these_rows <- .rowsToDrop(cg_table, gc_table, vr_table)

        if (!is.null(drop_these_rows)) {
            cg_table <- cg_table[!(rownames(cg_table) %in% drop_these_rows), , drop = FALSE]
            gc_table <- gc_table[!(rownames(gc_table) %in% drop_these_rows), , drop = FALSE]
            vr_table <- vr_table[!(rownames(vr_table) %in% drop_these_rows), , drop = FALSE]
        }
    }

    return (list(cg_table=cg_table,
                 gc_table=gc_table,
                 vr_table=vr_table))
}

# helper
# assumes inputs have the same sorted order
.rowsToDrop <- function(cg, gc, vr) {
    if (nrow(cg) != nrow(gc) && nrow(cg) != nrow(vr)) {
        stop("Error: Number of reads is unequal across CG, GC, and VARIANT tables")
    }

    drop <- c()
    for (row in nrow(cg)) {
        b_cg <- sum(is.na(cg[row,])) == ncol(cg)
        b_gc <- sum(is.na(gc[row,])) == ncol(gc)
        b_vr <- sum(is.na(vr[row,])) == ncol(vr)

        if (b_cg && b_gc && b_vr) {
            drop <- c(drop, row)
        }
    }

    drop
}

# helper
.tabulateRLE <- function(gr, col, include_snps, type=c("cg", "gc", "variant")) {
    type = match.arg(type)

    # there can be duplicate read names if not collapsed to fragment
    # this can occur if reads 1 and 2 originate from the "same" strand
    # make read names unique ahead of time...
    gr$readname <- make.unique(gr$readname)

    keep <- switch(
        type,
        "cg" = c("M", "U"),
        "gc" = c("O", "S"),
        "variant" = c("A", "T", "G", "C", "R", "Y")
    )
    if (include_snps && type != "variant") {
        keep <- c(keep, c("A", "T", "G", "C"))
    }

    # iterate through each read and decompose to position
    readlvl_gr <- do.call(
        "c",
        lapply(
            X = 1:length(gr),
            FUN = function(x) {
                sub_gr <- gr[x]

                # generate a per base array
                rle_vec <- unlist(strsplit(mcols(sub_gr)[, col], split = ""))
                names(rle_vec) <- seq(start(sub_gr), end(sub_gr))

                # filter out insertions
                rle_vec <- .filterInsertions(rle_vec)

                # keep values we want
                rle_vec_c <- rle_vec[rle_vec %in% keep]

                if (!length(rle_vec_c)) {
                    rle_c_df <- data.frame(chr = seqnames(sub_gr),
                                           start = start(sub_gr),
                                           end = start(sub_gr),
                                           strand = "*",
                                           base_status = NA,
                                           read_id = sub_gr$readname)
                    return(makeGRangesFromDataFrame(rle_c_df, keep.extra.columns = TRUE))
                }

                # turn back into GRanges
                rle_c_df <- data.frame(chr = seqnames(sub_gr),
                                       start = names(rle_vec_c),
                                       end = names(rle_vec_c),
                                       strand = sub_gr$bsstrand,
                                       base_status = rle_vec_c,
                                       read_id = sub_gr$readname)
                return(makeGRangesFromDataFrame(rle_c_df, keep.extra.columns = TRUE))
            }
        )
    )

    # find the dimension of collapsed values to fill in the matrix
    readlvl_gr_len <- length(unique(readlvl_gr))

    # make an empty matrix to fill in
    readlvl_emp_mat <- matrix(data = NA, nrow = length(unique(readlvl_gr$read_id)), ncol = readlvl_gr_len)
    rownames(readlvl_emp_mat) <- unique(readlvl_gr$read_id)
    colnames(readlvl_emp_mat) <- as.character(granges(unique(readlvl_gr)))

    # go by read and extract out methylation states
    readlvl_gr_mat <- as.matrix(cbind(as.character(granges(readlvl_gr)), readlvl_gr$read_id, readlvl_gr$base_status))
    readlvl_emp_mat[readlvl_gr_mat[,c(2,1)]] <- readlvl_gr_mat[,3]

    readlvl_emp_mat <- readlvl_emp_mat[,colSums(is.na(readlvl_emp_mat))<nrow(readlvl_emp_mat)]

    return(readlvl_emp_mat)
}

# helper to filter out insertions
# this is needed to reset coordinates properly
.filterInsertions <- function(readlvl_vec) {
    # the input here is a named vec of positions
    # we need to pull everything out that is not a lower case a,c,g,t,i
    # we can correct for new starts if someone has not filtered the first few bases
    exclude_bases <- c("a", "c", "g", "t", "i")

    # note: if a SNP has a base, it will be upper case
    # case sensitivity matters here
    # grab the start from the original string
    strt <- names(readlvl_vec)[1]
    filtrd_vec <- readlvl_vec[!readlvl_vec %in% exclude_bases]

    # short circuit if nothing is filtered
    if (suppressWarnings(all(names(filtrd_vec) == names(readlvl_vec), na.rm=TRUE))) {
        return(readlvl_vec)
    }

    # if the first base is no longer equal to the start from original read, reset to new start
    if (strt != names(filtrd_vec)[1]) {
        strt <- names(filtrd_vec)[1]
    }
    names(filtrd_vec) <- seq(as.numeric(strt), c(as.numeric(strt)+length(filtrd_vec)-1))

    return(filtrd_vec)
}

# helper to only keep sites within a given region of interest
.filterToRegion <- function(mat,
                            region = NULL) {

    # return coord sorted matrix if region == NULL
    if (is.null(region)) return(mat[,str_sort(colnames(mat), numeric=TRUE)])

    if (!is(region, "GRanges")) {
        # attempt to parse the standard chr#:start-end
        if (!grepl("\\:", region) & grepl("\\-", region)) {
            message("Not sure how to parse ", region)
            message("region should either be a GRanges of a specfic region or")
            stop("region should look something like 'chr6:1555-1900'")
        }
        chr <- strsplit(region, ":")[[1]][1]
        coords <- strsplit(region, ":")[[1]][2]
        strt <- strsplit(coords, "-")[[1]][1]
        end <- strsplit(coords, "-")[[1]][2]
        pos_to_include <- paste0(chr, ":", seq(strt, end))
    } else {
        # it's possible that multiple regions could be supplied...
        if (length(region > 1) & is(region, "GRanges")) {
            pos_to_include <- do.call(
                c,
                lapply(
                    1:length(region),
                    function(r) {
                        region.sub <- region[r]
                        return(paste0(seqnames(region.sub), ":", seq(start(region), end(region))))
                    }
                )
            )
        } else {
            # this is if a single region is supplied as a GRanges
            stopifnot(is(region, "GRanges"))
            pos_to_include <- paste0(seqnames(chr), ":", seq(start(region), end(region)))
        }
    }

    # subset to positions to include
    mat.sub <- mat[,substr(colnames(mat), 1, nchar(colnames(mat))-2) %in% pos_to_include]

    # order by chromosome position
    mat.sub <- mat.sub[,str_sort(colnames(mat.sub), numeric=TRUE)]

    return(mat.sub)
}

#' Plot the results of tabulateEpibed() as a quasi-lollipop plot.
#'
#' NOTE: If you run tabulateEpibed() in NOMe mode with include_empty_reads = TRUE,
#' then it is highly suggested you run plotEpiread() with show_all_points = TRUE also.
#' If not, it is possible that reads in the CG and GC plots will not line up, depending
#' on if there are empty reads in one plot, but not the other.
#'
#' @param mat Input matrix that comes out of tabulateEpibed()
#' @param plot_read_avg Whether to also plot the average methylation state (default: TRUE)
#' @param show_readnames Whether to show the read names (default: TRUE)
#' @param show_positions Whether to show the genomic positions (default: TRUE)
#' @param show_all_points Whether to show all points (i.e., empty reads, unknown, and filtered sites) (default: FALSE)
#' @param meth_color What color should the methylated states be (default: '#FF2400')
#' @param unmeth_color What color should the unmethylated states be (default: '#6495ED')
#' @param na_color What color should the NA values be (default: 'grey')
#' @param size What size the points should be drawn as (default: 6)
#'
#' @return An epiread ggplot object or list of ggplot objects if plot_read_avg is TRUE
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#'
#' @export
#'
#' @examples
#'
#' epibed.nome <- system.file("extdata", "hct116.nome.epiread.gz",
#'                            package="biscuiteer")
#' epibed.nome.gr <- readEpibed(epibed = epibed.nome, is.nome = TRUE,
#'                              genome = "hg19", chr = "chr1")
#' epibed.tab.nome <- tabulateEpibed(epibed.nome.gr)
#' plotEpiread(epibed.tab.nome$gc_table)
#'
plotEpiread <- function(mat,
                        plot_read_avg = TRUE,
                        show_readnames = TRUE,
                        show_positions = TRUE,
                        show_all_points = FALSE,
                        meth_color = "#FF2400",
                        unmeth_color = "#6495ED",
                        na_color = "grey",
                        size = 6) {

    # check if input is a matrix
    if (is.list(mat) & is(mat[[1]], "matrix")) {
        message("Input given is likely the output of tabulateEpibed() in NOMe mode.")
        stop("Please select either $cg_table, $gc_table, or $vr_table to plot.")
    }
    if (!is.list(mat) & !is(mat, "matrix")) {
        stop("Input needs to be a matrix generated by tabulateEpibed()")
    }

    # If running in NOMe mode and include_empty_reads = TRUE in tabulateEpibed, then the reads in the CG and GC plots
    # may not line up. To ensure alignment, it's best to set show_all_points = TRUE in this instance.
    mat.melt <- .makePlotData(mat, show_all_points)

    # plot epiread
    ql_theme <- .set_ql_theme(show_readnames, show_positions)
    plt <- .epiClustPlot(mat.melt, ql_theme, meth_color, unmeth_color, na_color, size)

    # average methylation
    # FIXME: If having SNP plot separate, find way to avoid making average plot for SNPs
    if (plot_read_avg) {
        plt_avg <- .plotAvg(mat, meth_color, unmeth_color, ql_theme, size)
        return(list(epistate=plt,
                    meth_avg=plt_avg))
    } else {
        return(plt)
    }
}

# helper to make the melted dataset for plotting
.makePlotData <- function(mat, show_all_points) {
    # auto-detect input type
    is.cg = FALSE
    is.gc = FALSE
    is.vr = FALSE
    if(any(mat %in% c("M", "U"))) is.cg = TRUE
    if(any(mat %in% c("S", "O"))) is.gc = TRUE
    if(any(mat %in% c("A", "C", "T", "G", "R", "Y"))) is.vr = TRUE

    # break if both are FALSE
    if (isFALSE(is.cg) & isFALSE(is.gc) & isFALSE(is.vr)) {
        message("Don't know what to do with input.")
        stop("Please run tabulateEpibed() first to produce input for this function.")
    }

    # cast to a 'melted' data frame
    mat.melt <- reshape2::melt(mat, id.vars = rownames(mat))
    if (!show_all_points) {
        mat.melt <- subset(mat.melt, !is.na(value))
    }

    return(mat.melt)
}

# helper to calculate avg methylation of a region and plot it
.plotAvg <- function(mat, meth_color, unmeth_color, theme, size) {
    # check if input is a matrix
    if (!is(mat, "matrix")) {
        stop("Input needs to be a matrix")
    }

    mat[mat %in% c("M", "O")] <- 1
    mat[mat %in% c("U", "S")] <- 0
    mat[mat %in% c("A", "T", "G", "C")] <- -1

    mat <- apply(mat, 2, as.numeric)
    mat.meth.avg <- data.frame(avg_meth = colMeans(mat, na.rm = TRUE))
    mat.meth.avg <- na.omit(mat.meth.avg) # remove empty rows that are somehow making it all the way here
    mat.meth.avg$avg_meth <- ifelse(mat.meth.avg$avg_meth < 0, NA, mat.meth.avg$avg_meth)

    mat.meth.avg$position <- rownames(mat.meth.avg)
    mat.meth.avg$y <- ifelse(is.na(mat.meth.avg$avg_meth), "SNP status", "Average methylation status")

    plt_avg <- ggplot(mat.meth.avg, aes(x = position, y = y)) +
        geom_point(aes(fill = avg_meth), size=size, pch=21, color="black", na.rm=TRUE) +
        scale_fill_gradient(low = unmeth_color, high = meth_color, limits = c(0,1)) +
        scale_x_discrete(name ="", limits=rownames(mat.meth.avg)) +
        scale_y_discrete(limits = c("Average methylation status")) +
        guides(color = "legend") +
        theme

    return(plt_avg)
}

# helper to set the ql theme depending on readnames and position toggles
.set_ql_theme <- function(show_readnames = TRUE, show_positions = TRUE) {
    ql_theme <- theme_bw(12) +
        theme(
            axis.title = element_blank(),
            legend.title = element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(color = "black")
        )

    if (!show_readnames) {
        # set the theme
        ql_theme <- ql_theme +
            theme(
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank()
            )
    }

    if (!show_positions) {
        ql_theme <- ql_theme +
            theme(
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank()
            )
    }

    return(ql_theme)
}

# helper to plot epireads given a matrix of read, pos, value
# TODO: Remove default options as this is a non-exported function
#       Will need to then define these values for the .epiClustPlot call in the epistateCaller function
.epiClustPlot <- function(matplotdata,
                          ql_theme,
                          meth_color = "#FF2400",
                          unmeth_color = "#6495ED",
                          na_color = "grey",
                          size = 6) {

    snp_list <- c("A", "T", "G", "C", "R", "Y")

    plt <- ggplot(matplotdata, aes(x = Var2, y = Var1)) +
        geom_point(aes(fill = value), size=size, pch=21, color="black") +
        guides(color = "legend") +
        geom_label(data = subset(matplotdata, value %in% snp_list), aes(label = value)) +
        ql_theme

    if (any(matplotdata$value %in% c("M", "U"))) {
        plt <- plt +
            scale_fill_manual(values = c(M = meth_color, U = unmeth_color), na.value = na_color)
    } else {
        plt <- plt +
            scale_fill_manual(values = c(O = meth_color, S = unmeth_color), na.value = na_color)
    }

    return(plt)
}

#' Calculate and plot a hierarchical clustering of the epireads
#'
#' WARNING: This is a beta function and should not be used for publication level analyses!!!!
#'
#' @param mat Input matrix that comes out of tabulateEpibed()
#' @param stringdist_method stringdist::stringdist algorithm to use (default: "hamming")
#' @param hclust_method Clustering algorithm to use (default: "ward.D2")
#' @param plot Whether to plot the clustered epireads (default: TRUE)
#'
#' @return A hierarchical cluster (hclust) object
#'
#' @import ggtree
#' @importFrom stringdist stringdist
#' @importFrom cowplot plot_grid
#'
#' @examples
#'
#' #epibed.nome <- system.file("extdata", "hct116.nome.epiread.gz",
#' #                           package="biscuiteer")
#' #epibed.nome.gr <- readEpibed(epibed = epibed.nome, is.nome = TRUE,
#' #                             genome = "hg19", chr = "chr1")
#' #epibed.tab.nome <- tabulateEpibed(epibed.nome.gr)
#' #epistateCaller(epibed.tab.nome)
#'
epistateCaller <- function(mat,
                           stringdist_method="hamming",
                           hclust_method="ward.D2",
                           plot = TRUE) {
    # check for both cg and gc tables
    if (is.list(mat)) {
        stopifnot(exists("cg_table", where = mat) & exists("gc_table", where = mat))

        mat.merge <- merge(mat$cg_table, mat$gc_table, by = 0)
        row.names(mat.merge) <- mat.merge[, 1]
        mat.merge <- mat.merge[-1]
    } else {
        mat.merge <- mat
    }

    mat.merge[is.na(mat.merge) |
              mat.merge == "A" | mat.merge == "T" |
              mat.merge == "G" | mat.merge == "C" |
              mat.merge == "U" | mat.merge == "S"] <- 0
    mat.merge[mat.merge == "M" | mat.merge == "O"] <- 1

    epistates <- apply(mat.merge, 1, paste0, collapse = "")
    hamming <- outer(epistates, epistates, stringdist, method = stringdist_method)

    hcluster <- hclust(as.dist(hamming), method = hclust_method)
    if (!plot) {
        return(hcluster)
    }

    if (!is.list(mat)) {
        matplotdata <- list(.makePlotData(mat, show_all_points = FALSE))
    } else {
        matplotdata <- lapply(
            mat, function(meth_table) { .makePlotData(meth_table, show_all_points = FALSE) }
        )
    }

    tree <- ggtree(as.dendrogram(hcluster), branch.length = "none")

    ql_theme <- .set_ql_theme(TRUE, TRUE)

    plots <- lapply(
        matplotdata,
        function(m) {
            m$Var1 <- factor(m$Var1, levels = rev(get_taxa_name(tree)))
            plt <- .epiClustPlot(m, ql_theme)

            clust_plot <- plot_grid(
                tree, NULL, plt,
                rel_widths = c(1, -0.05, 2), nrow = 1, align = 'h'
            )
            return(clust_plot)
        }
    )

    return(plots)
}
