#' Tabulate to a read or fragment level matrix for CGH or GCH context methylation.
#'
#' The input for this function is the GRanges output from biscuiteer::readEpibed().
#' It will auto-detect whether this a BS-seq or a NOMe-seq run and tabulate accordingly.
#'
#' @param gr The epibed GRanges object from readEpibed()
#' @param region Either a GRanges of regions to subset to or explicit region (e.g., chr6:1555-1900)
#' @param include_empty_reads Whether to include reads that contain no methylated (or SNP) sites (default: FALSE)
#' @param include_snps Whether to include SNPs or not, only matters for epiBEDs from BISCUIT v1.0 or 1.1 (default: TRUE)
#' @param merge_cg_vr Whether to create a table with the CpG and SNP tables merged (default: TRUE)
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
#' epibed.nome <- system.file("extdata", "hct116.nome.epibed.gz",
#'                            package="biscuiteer")
#' epibed.nome.gr <- readEpibed(epibed = epibed.nome,
#'                              genome = "hg19", chr = "chr1")
#' epibed.tab.nome <- tabulateEpibed(epibed.nome.gr)
#'
tabulateEpibed <- function(gr,
                           region = NULL,
                           include_empty_reads = FALSE,
                           include_snps = TRUE,
                           merge_cg_vr = TRUE) {

    # check if a GRanges
    stopifnot(is(gr, "GRanges"))

    # tabulate strings
    #const char SKIP_EPI = '-';
    #const char SKIP_INS = 'i';
    #const char SKIP_DEL = 'd';
    #const char FILTERED = 'F';
    #const char IGNORED  = 'x';
    #const char DELETION = 'D';
    #const char SOFTCLIP = 'P';
    #const char METHYLAT = 'M';
    #const char UNMETHYL = 'U';
    #const char OPEN_ACC = 'O';
    #const char SHUT_ACC = 'S';
    #const char AMBIG_GA = 'R';
    #const char AMBIG_CT = 'Y';

    # we need to build up a table of CGH and GCH
    if (!all(is.na(gr$CG_decode))) {
        cg_table <- .tabulateRLE(gr, col = "CG_decode", include_snps = include_snps, type = "cg")
        cg_table <- .filterToRegion(cg_table, region = region)
        cg_table <- .mergeColumns(cg_table, type="cg")
    } else {
        cg_table <- matrix(data = NA, nrow = length(gr$readname), ncol = 1)
        rownames(cg_table) <- gr$readname
        colnames(cg_table) <- "cg_missing"
    }

    if (!all(is.na(gr$GC_decode))) {
        gc_table <- .tabulateRLE(gr, col = "GC_decode", include_snps = include_snps, type = "gc")
        gc_table <- .filterToRegion(gc_table, region = region)
        gc_table <- .mergeColumns(gc_table, type="gc")
    } else {
        gc_table <- matrix(data = NA, nrow = length(gr$readname), ncol = 1)
        rownames(gc_table) <- gr$readname
        colnames(gc_table) <- "gc_missing"
    }

    if (!all(is.na(gr$VAR_decode))) {
        vr_table <- .tabulateRLE(gr, col = "VAR_decode", include_snps = include_snps, type = "variant")
        vr_table <- .filterToRegion(vr_table, region = region)
        vr_table <- .mergeColumns(vr_table, type="snp")
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

    if (merge_cg_vr) {
        tab_cgvr = .combineCgVr(cg_table, vr_table)
    } else {
        tab_cgvr = NULL
    }

    return (list(cg_table=cg_table,
                 gc_table=gc_table,
                 vr_table=vr_table,
                 tab_cgvr=tab_cgvr))
}

# helper
# fname - name of column whose entries will be first in paste
# lname - name of column whose entries will be last in paste
.mergeWithNAs <- function(mat, fname, lname) {
    tmp <- paste0(
        replace(mat[, fname], is.na(mat[, fname]), "N"),
        replace(mat[, lname], is.na(mat[, lname]), "N")
    )

    tmp
}

# helper
.combineCgVr <- function(cg, vr) {
    comb <- cbind(cg, vr)
    comb <- comb[, sort(colnames(comb))]

    to_drop <- c()
    for (i in seq_along(colnames(comb))) {
        creg = .parseRegion(colnames(comb)[i])
        cnam = colnames(comb)[i]

        if (creg[["end"]] - creg[["start"]] > 0) {
            if (i == 1) {
                # only need to check next column
                nreg = .parseRegion(colnames(comb)[i+1])
                nnam = colnames(comb)[i+1]

                if ((nreg[["end"]] - nreg[["start"]] == 0) && (creg[["end"]] == nreg[["start"]])) {
                    tmp <- .mergeWithNAs(comb, cnam, nnam)
                    comb[, cnam] <- replace(tmp, tmp == "NN", NA)
                    comb[, nnam] <- NA
                    to_drop <- c(to_drop, nnam)
                } else {
                    next
                }
            } else if (i == ncol(comb)) {
                # only to check previous column
                preg = .parseRegion(colnames(comb)[i-1])
                pnam = colnames(comb)[i-1]

                if ((preg[["end"]] - preg[["start"]] == 0) && (preg[["end"]] == creg[["start"]])) {
                    tmp <- .mergeWithNAs(comb, pnam, cnam)
                    comb[, cnam] <- replace(tmp, tmp == "NN", NA)
                    comb[, pnam] <- NA
                    to_drop <- c(to_drop, pnam)
                } else {
                    next
                }
            } else {
                # need to check both previous and last columns
                preg = .parseRegion(colnames(comb)[i-1])
                pnam = colnames(comb)[i-1]

                if ((preg[["end"]] - preg[["start"]] == 0) && (preg[["end"]] == creg[["start"]])) {
                    tmp <- .mergeWithNAs(comb, pnam, cnam)
                    comb[, cnam] <- replace(tmp, tmp == "NN", NA)
                    comb[, pnam] <- NA
                    to_drop <- c(to_drop, pnam)
                }

                nreg = .parseRegion(colnames(comb)[i+1])
                nnam = colnames(comb)[i+1]

                if ((nreg[["end"]] - nreg[["start"]] == 0) && (creg[["end"]] == nreg[["start"]])) {
                    tmp <- .mergeWithNAs(comb, cnam, nnam)
                    comb[, cnam] <- replace(tmp, tmp == "NN", NA)
                    comb[, nnam] <- NA
                    to_drop <- c(to_drop, nnam)
                }
            }
        } else {
            next
        }
    }

    comb <- comb[, !(colnames(comb) %in% to_drop)]

    comb
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
        par <- .parseRegion(region)
        pos_to_include <- paste0(par[["chr"]], ":", seq(par[["start"]], par[["end"]]))
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
            pos_to_include <- paste0(seqnames(region), ":", seq(start(region), end(region)))
        }
    }

    # subset to positions to include
    mat.sub <- mat[,substr(colnames(mat), 1, nchar(colnames(mat))-2) %in% pos_to_include]

    # order by chromosome position
    mat.sub <- mat.sub[,str_sort(colnames(mat.sub), numeric=TRUE)]

    return(mat.sub)
}

# helper
.parseRegion <- function(region) {
    chr    <- strsplit(region, ":")[[1]][1]
    coords <- strsplit(region, ":")[[1]][2]
    start  <- strsplit(coords, "-")[[1]][1]
    end    <- strsplit(coords, "-")[[1]][2]

    return(list(chr=chr, start=as.integer(start), end=as.integer(end)))
}

# helper
# TODO: Refactor this behomoth
.mergeColumns <- function(mat, type=c("cg", "gc", "snp")) {
    type <- match.arg(type)

    prev <- NULL
    to_drop <- c()
    if (type == "cg") {
        for (i in seq_len(length(colnames(mat)))) {
            cname <- colnames(mat)[i]
            if (cname == "cg_missing") {
                return(mat)
            }
            len <- nchar(cname)
            if (is.null(prev)) {
                last_char <- substr(cname, len, len)
                if (last_char == "+") {
                    prev <- cname
                    next
                } else if (last_char == "*") {
                    cat("* in a CG found\n")
                    break
                } else {
                    pieces <- .parseName(cname)
                    new_cname <- paste(pieces$chr, pieces$loc-1, "+", sep=":")
                    colnames(mat)[i] <- new_cname
                    next
                }
            }

            ppieces <- .parseName(prev)
            cpieces <- .parseName(cname)

            if (ppieces$loc+1 == cpieces$loc) {
                for (rname in rownames(mat)) {
                    if (!is.na(mat[rname, prev]) && !is.na(mat[rname, cname])) {
                        message("Malformed file: Found values defined in both OT and OB strands!")
                        stop("Read name: ", rname, ", location: ", prev, ", ", cname)
                    } else if (is.na(mat[rname, prev]) && !is.na(mat[rname, cname])) {
                        mat[rname, prev] <- mat[rname, cname]
                        mat[rname, cname] <- NA
                    } else {
                        # Either both are NAs or "-" strand is, but not "+" strand
                        next
                    }
                }

                if (all(is.na(mat[cname]))) {
                    to_drop <- c(to_drop, cname)
                }
            } else {
                last_char <- substr(cname, len, len)
                if (last_char == "+") {
                    prev <- cname
                    next
                } else {
                    new_cname <- paste(cpieces$chr, cpieces$loc-1, "+", sep=":")
                    colnames(mat)[i] <- new_cname
                }
            }

            prev <- NULL
        }

        mat <- mat[,!(colnames(mat) %in% to_drop)]

        for (i in seq_len(length(colnames(mat)))) {
            pieces <- .parseName(colnames(mat)[i])
            colnames(mat)[i] <- paste0(pieces$chr, ":", pieces$loc, "-", pieces$loc+1)
        }
    } else if (type == "gc") {
        for (i in seq_len(length(colnames(mat)))) {
            cname <- colnames(mat)[i]
            if (cname == "gc_missing") {
                return(mat)
            }
            len <- nchar(cname)
            if (is.null(prev)) {
                last_char <- substr(cname, len, len)
                if (last_char == "-") {
                    prev <- cname
                    next
                } else if (last_char == "*") {
                    cat("* in a GC found\n")
                    break
                } else {
                    pieces <- .parseName(cname)
                    new_cname <- paste(pieces$chr, pieces$loc-1, "-", sep=":")
                    colnames(mat)[i] <- new_cname
                    next
                }
            }

            ppieces <- .parseName(prev)
            cpieces <- .parseName(cname)

            if (ppieces$loc+1 == cpieces$loc) {
                for (rname in rownames(mat)) {
                    if (!is.na(mat[rname, prev]) && !is.na(mat[rname, cname])) {
                        message("Malformed file: Found values defined in both OT and OB strands!")
                        stop("Read name: ", rname, ", location: ", prev, ", ", cname)
                    } else if (is.na(mat[rname, prev]) && !is.na(mat[rname, cname])) {
                        mat[rname, prev] <- mat[rname, cname]
                        mat[rname, cname] <- NA
                    } else {
                        # Either both are NAs or "-" strand is, but not "+" strand
                        next
                    }
                }

                if (all(is.na(mat[cname]))) {
                    to_drop <- c(to_drop, cname)
                }
            } else {
                last_char <- substr(cname, len, len)
                if (last_char == "-") {
                    prev <- cname
                    next
                } else {
                    new_cname <- paste(cpieces$chr, cpieces$loc-1, "-", sep=":")
                    colnames(mat)[i] <- new_cname
                }
            }

            prev <- NULL
        }

        mat <- mat[,!(colnames(mat) %in% to_drop)]

        for (i in seq_len(length(colnames(mat)))) {
            pieces <- .parseName(colnames(mat)[i])
            colnames(mat)[i] <- paste0(pieces$chr, ":", pieces$loc, "-", pieces$loc+1)
        }
    } else {
        for (i in seq_len(length(colnames(mat)))) {
            cname <- colnames(mat)[i]
            if (cname == "vr_missing") {
                return(mat)
            }
            len <- nchar(cname)
            if (is.null(prev)) {
                last_char <- substr(cname, len, len)
                if (last_char == "-") {
                    prev <- cname
                    next
                } else if (last_char == "*") {
                    cat("* in a SNP found\n")
                    break
                } else {
                    pieces <- .parseName(cname)
                    new_cname <- paste(pieces$chr, pieces$loc, "-", sep=":")
                    colnames(mat)[i] <- new_cname
                    next
                }
            }

            ppieces <- .parseName(prev)
            cpieces <- .parseName(cname)

            if (ppieces$loc == cpieces$loc) {
                for (rname in rownames(mat)) {
                    if (!is.na(mat[rname, prev]) && !is.na(mat[rname, cname])) {
                        message("Malformed file: Found values defined in both OT and OB strands!")
                        stop("Read name: ", rname, ", location: ", prev, ", ", cname)
                    } else if (is.na(mat[rname, prev]) && !is.na(mat[rname, cname])) {
                        mat[rname, prev] <- mat[rname, cname]
                        mat[rname, cname] <- NA
                    } else {
                        # Either both are NAs or "-" strand is, but not "+" strand
                        next
                    }
                }

                if (all(is.na(mat[cname]))) {
                    to_drop <- c(to_drop, cname)
                }
            } else {
                last_char <- substr(cname, len, len)
                if (last_char == "-") {
                    prev <- cname
                    next
                } else {
                    new_cname <- paste(cpieces$chr, cpieces$loc, "-", sep=":")
                    colnames(mat)[i] <- new_cname
                }
            }

            prev <- NULL
        }

        mat <- mat[,!(colnames(mat) %in% to_drop)]

        for (i in seq_len(length(colnames(mat)))) {
            pieces <- .parseName(colnames(mat)[i])
            colnames(mat)[i] <- paste0(pieces$chr, ":", pieces$loc, "-", pieces$loc)
        }
    }

    mat
}

# helper
.parseName <- function(cname) {
    pieces <- strsplit(cname, ":")[[1]]

    list(
        chr = pieces[1],
        loc = as.integer(pieces[2]),
        stn = pieces[3]
    )
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

