#' Tabulate to a read or fragment level matrix for
#' CGH or GCH context methylation. The input for this
#' function is the GRanges output from readEpibed().
#' It will auto-detect whether this a BS-seq or a
#' NOMe-seq run and tabulate accordingly.
#'
#' @param gr The epibed GRanges object from readEpibed()
#' @param region Either a GRanges of regions to subset to or explicit region (e.g. chr6:1555-1900)
#' @param filter_empty_reads Whether to filter out reads that contain no methylated sites (default is TRUE)
#' 
#' @return A matrix or list of matrices
#' @export
#' 
#' @import GenomicRanges
#'
#' @examples
#' epibed.nome <- system.file("extdata", "hct116.nome.epiread",
#'                            package="biscuiteer")
#' epibed.nome.gr <- readEpibed(epibed = epibed.nome, is.nome = TRUE,
#'                              genome = "hg19")
#' epibed.tab.nome <- tabulateEpibed(epibed.nome.gr)

tabulateEpibed <- function(gr,
                           region = NULL,
                           filter_empty_reads = TRUE) {
  
  # check if a GRanges
  stopifnot(is(gr, "GRanges"))
  
  # autodetect if it's a NOMe epibed
  is.nome = FALSE
  if ("GC_decode" %in% names(mcols(gr))) is.nome = TRUE
  
  # tabulate strings
  # char filtered = 'F';
  # char ignored  = 'x';
  # char deletion = 'D';
  # char softclip = 'P';
  # char methylat = 'M';
  # char unmethyl = 'U';
  # char open_acc = 'O';
  # char shut_acc = 'S';
  
  # we need to build up a table of CGH and GCH
  # cg_table <- .tabulateCGH(gr)
  cg_table <- .tabulateRLE(gr, cg = TRUE)
  cg_table <- .filterToRegion(cg_table,
                            region = region)
  if (filter_empty_reads) {
    cg_table <- .filterEmptyReads(cg_table)
  }
  if (is.nome) {
    gc_table <- .tabulateRLE(gr, cg = FALSE)
    gc_table <- .filterToRegion(gc_table,
                              region = region)
    if (filter_empty_reads) {
      gc_table <- .filterEmptyReads(gc_table)
    }
    return(list(cg_table = cg_table,
                gc_table = gc_table))
  }
  
  return(cg_table)
}

# helper to tabulate to a CGH or GCH table
.tabulateRLE <- function(gr, cg = TRUE) {
  # there can be duplicate read names if not collapsed to fragment
  # this can occur if reads 1 and 2 originate from the "same" strand
  # make read names unique ahead of time...
  gr$readname <- make.unique(gr$readname)
  
  # iterate through each read and decompose to position
  readlvl_gr <- do.call("c", lapply(X = 1:length(gr),
                                    FUN = function(x) {
    sub_gr <- gr[x]
    # generate a per base array
    pos_vec <- seq(start(sub_gr), end(sub_gr))
    if (cg) {
      rle_vec <- unlist(strsplit(sub_gr$CG_decode, split = ""))
    } else {
      rle_vec <- unlist(strsplit(sub_gr$GC_decode, split = ""))
    }
    names(rle_vec) <- pos_vec
    # filter out insertions
    rle_vec <- .filterInsertions(rle_vec)
    # keep the C status
    if (cg) {
      rle_vec_c <- rle_vec[rle_vec %in% c("M", "U")]
    } else {
      rle_vec_c <- rle_vec[rle_vec %in% c("O", "S")]
    }
    if (!length(rle_vec_c)) {
      return()
    }
    # turn back into GRanges
    rle_c_df <- data.frame(chr = seqnames(sub_gr),
                            start = names(rle_vec_c),
                            end = names(rle_vec_c),
                            meth_status = rle_vec_c,
                            read_id = sub_gr$readname)
    return(makeGRangesFromDataFrame(rle_c_df,
                                    keep.extra.columns = TRUE))
  }))
  
  # find the dimension of collapsed Cs to fill in the matrix
  readlvl_gr_len <- length(unique(readlvl_gr))
  # make an empty matrix to fill in 
  readlvl_emp_mat <- matrix(data = NA,
                            nrow = length(unique(readlvl_gr$read_id)),
                            ncol = readlvl_gr_len)
  rownames(readlvl_emp_mat) <- unique(readlvl_gr$read_id)
  colnames(readlvl_emp_mat) <- as.character(granges(unique(readlvl_gr)))
  
  # go by read and extract out methylation states
  readlvl_gr_mat <- as.matrix(cbind(as.character(granges(readlvl_gr)),
                                    readlvl_gr$read_id,
                                    readlvl_gr$meth_status))
  readlvl_emp_mat[readlvl_gr_mat[,c(2,1)]] <- readlvl_gr_mat[,3]
  return(readlvl_emp_mat)
}

# helper to filter out indels and softclips
# this is needed to reset coordinates properly
.filterInsertions <- function(readlvl_vec) {
  # the input here is a named vec of positions
  # we need to pull everything out that is not
  # a lower case a,c,g,t
  # we can correct for new starts if someone has
  # not filtered the first few bases
  exclude_bases <- c("a", "c",
                     "g", "t")
  # note: if a SNP has a base, it will be upper case
  # case sensitivity matters here
  # grab the start from the original string
  strt <- names(readlvl_vec)[1]
  filtrd_vec <- readlvl_vec[!readlvl_vec %in% exclude_bases]
  # short circuit if nothing is filtered
  if (suppressWarnings(all(names(filtrd_vec) == names(readlvl_vec)))) {
    return(readlvl_vec)
  }
  # if the first base is no longer equal to the start
  # from original read, reset to new start
  if (strt != names(filtrd_vec)[1]) {
    strt <- names(filtrd_vec)[1]
    }
  names(filtrd_vec) <- seq(as.numeric(strt),
                           c(as.numeric(strt)+length(filtrd_vec)-1))
  return(filtrd_vec)
}

# helper to only keep sites within a given region of interest
.filterToRegion <- function(mat,
                            region = NULL) {
  if (is.null(region)) return(mat)
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
    pos_to_include <- paste0(chr, ":",
                             seq(strt, end))
  } else {
    # it's possible that multiple regions could be supplied...
    if (length(region > 1) & is(region, "GRanges")) {
      pos_to_include <- do.call(c, lapply(1:length(region), function(r) {
        region.sub <- region[r]
        return(paste0(seqnames(region.sub), ":",
                      seq(start(region),
                          end(region))))
      }))
    } else {
      # this is if a single region is supplied as a GRanges
      stopifnot(is(region, "GRanges"))
      pos_to_include <- paste0(seqnames(chr), ":",
                               seq(start(region),
                                   end(region)))
    }
  }
  # subset
  mat.sub <- mat[,colnames(mat) %in% pos_to_include]
  # order
  mat.sub <- mat.sub[,order(colnames(mat.sub))]
  return(mat.sub)
}

# helper to remove empty reads
# an 'empty' read is one with all NAs
.filterEmptyReads <- function(mat) {
  # input is a matrix after tabulateEpiread is done
  # filter reads
  mat.sub <- mat[rowMeans(is.na(mat)) < 1,]
  # order
  mat.sub <- mat.sub[,order(colnames(mat.sub))]
  return(mat.sub)
}

#' Plot the results of tabulateEpibed() as a 
#' quasi-lollipop plot.
#'
#' @param mat Input matrix that comes out of tabulateEpibed()
#' @param plot_read_ave Whether to also plot the average methylation state (default: TRUE)
#' @param show_readnames Whether to show the read names (default: TRUE)
#' @param show_positions Whether to show the genomic positions (default: TRUE)
#' @param meth_color What color should the methylated states be (default: 'black')
#' @param unmeth_color What color should the unmethylated states be (default: 'white')
#' @param na_color What color should the NA values be (default: 'darkgray')
#' @param background_color What color the background of the plot should be (default: '#A3D0E9')
#' 
#' @return An epiread ggplot object or list of ggplot objects if plot_read_ave is TRUE
#' 
#' @import ggplot2
#' @importFrom reshape2 melt
#' 
#' @export
#'
#' @examples
#' epibed.nome <- system.file("extdata", "hct116.nome.epiread",
#'                            package="biscuiteer")
#' epibed.nome.gr <- readEpibed(epibed = epibed.nome, is.nome = TRUE,
#'                              genome = "hg19")
#' epibed.tab.nome <- tabulateEpibed(epibed.nome.gr)
#' plotEpiread(epibed.tab.nome$gc_table)

plotEpiread <- function(mat, plot_read_ave = TRUE,
                        show_readnames = TRUE,
                        show_positions = TRUE,
                        unmeth_color = "white",
                        meth_color = "black",
                        na_color = "darkgray",
                        background_color = "#A3D0E9") {
  
  # check if input is a matrix
  if (!is(mat, "matrix")) {
    stop("Input needs to be a matrix generated by tabulateEpibed()")
  }
  
  # set themes
  ql_theme <- theme_bw(12) + theme(
    axis.title = element_blank(),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45,
                               vjust = 1,
                               hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "black"),
    panel.background = element_rect(fill = background_color))
  
  if (!show_readnames) {
    # set the theme
    ql_theme <- ql_theme + theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank())
  }
  
  if (!show_positions) {
    ql_theme <- ql_theme + theme(axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank())
  }
  
  # auto-detect input type
  is.cg = FALSE
  is.gc = FALSE
  if(any(mat %in% c("M", "U"))) is.cg = TRUE
  if(any(mat %in% c("S", "O"))) is.gc = TRUE
  
  # break if both are FALSE
  if (isFALSE(is.cg) & isFALSE(is.gc)) {
    message("Don't know what to do with input.")
    stop("Please run tabulateEpibed() first to produce input for this function.")
  }
  
  # compute average methylation state
  if (plot_read_ave) {
    # ick... copy in memory
    mat.binary <- mat
    mat.binary[mat.binary %in% c("M", "O")] <- 1
    mat.binary[mat.binary %in% c("U", "S")] <- 0
    mat.binary <- apply(mat.binary, 2, as.numeric)
    mat.meth.ave <- data.frame(ave_meth = colMeans(mat.binary, na.rm = TRUE))
    mat.meth.ave$position <- rownames(mat.meth.ave)
    mat.meth.ave$y <- "Average methylation status"
  }
  
  # cast to a 'melted' data frame
  mat.melt <- reshape2::melt(mat, id.vars = rownames(mat))
  
  # plot
  if (is.cg) {
    plt <- ggplot(mat.melt, aes(x = Var2, y = Var1)) +
      geom_point(aes(fill = value), size=6, pch=21, color="black") +
      scale_fill_manual(values = c(M=meth_color, U=unmeth_color),
                         na.value = na_color) +
      guides(color = "legend") +
      ql_theme
  } else {
    plt <- ggplot(mat.melt, aes(x = Var2, y = Var1)) +
      geom_point(aes(fill = value), size=6, pch=21, color="black") +
      scale_fill_manual(values = c(O=meth_color, S=unmeth_color),
                        na.value = na_color) +
      guides(color = "legend") +
      ql_theme
      
  }
  
  if (plot_read_ave) {
    plt_ave <- ggplot(mat.meth.ave, aes(x = position, y = y)) +
      geom_point(aes(fill = ave_meth), size=6, pch=21, color="black") +
      scale_fill_gradient(low = unmeth_color,
                          high = meth_color,
                          limits = c(0,1)) +
      guides(color = "legend") +
      ql_theme
    return(list(epistate=plt,
                meth_ave=plt_ave))
  } else {
    return(plt)
  }
}

