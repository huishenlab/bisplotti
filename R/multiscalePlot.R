#' Create the multiscale plot from a GRangesList
#'
#' @param grl GRangesList of a single set of multiscale BED files
#' @param colors RColorBrewer or viridis color scheme to use (default: YlGnBu)
#' @param na_color Color to give NA values (default: darkgray)
#' @param what Column to plot (default: score)
#' @param config path to the bin/config.yaml of the multiscale Snakemake
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
#' files.loc <- system.file("extdata", package="bisplotti")
#'
#' files <- lapply(
#'     list.files(files.loc, pattern="Heyn_2012_100yr", full.names=TRUE), function(x) {
#'         return(rtracklayer::import(x, format="bedGraph"))
#'     }
#' )
#' cnames <- list.files(files.loc, pattern="Heyn_2012_100yr")
#' cnames <- gsub("Heyn_2012_100yr_subset_", "", cnames)
#' cnames <- gsub(".bed.gz", "", cnames)
#' names(files) <- as.numeric(cnames)
#'
#' files.grl <- as(files, "GRangesList")
#'
#' multiscaleMethylationPlot(files.grl)
#'
multiscaleMethylationPlot <- function(grl,
                                      colors="YlGnBu",
                                      na_color="darkgray",
                                      what="score",
                                      config=NULL,
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
    if (!is.null(chr.to.plot)) {
        rownames(score.matrix) <- as.numeric(start(long.mat))
    } else {
        rownames(score.matrix) <- as.character(granges(long.mat))
    }

    score.matrix <- as.data.frame(t(score.matrix))
    score.matrix$res <- rownames(score.matrix)
    score.matrix.melt <- reshape2::melt(score.matrix, id.vars = "res")
    score.matrix.melt$res <- levels(factor(score.matrix.melt$res))

    # Create plot
    g <- ggplot(score.matrix.melt, aes(x = variable, y = res, fill = value)) +
        geom_raster(interpolate = TRUE) +
        theme_bw(12)

    # Add y-axis ticks (if config is provided)
    if (!is.null(config)) {
        steps <- .getMultiscaleSteps(config) # steps

        g <- g + 
            scale_y_discrete(breaks = unique(round(steps,0)), name="",
                             labels = paste((10^unique(round(steps,0)))/1e6,'Mb'))

    } else {
        g <- g +
        theme(
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
        )
    }

    # Add x-axis values if plotting a chromosome
    if (!is.null(chr.to.plot)) {
        vals <- as.numeric(as.character(levels(score.matrix.melt$variable)))
        g <- g + scale_x_discrete(breaks=seq(min(vals), max(vals), 1e7),
                                  labels=round(seq(min(vals)-1, max(vals)/1e6, 10))) +
            xlab("Position along Chromosome (Mb)")
    } else {
        g <- g +
            theme(
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
            )
    }

    # Add color to plot
    if (colors %in% c("viridis", "magma", "plasma", "inferno", "cividis", "mako", "rocket", "turbo")) {
        g <- g + scale_fill_viridis(option = colors, direction = 1, na.value = na_color,
                                    name = "Methylation Level", limits = c(0,1))
    } else {
        g <- g + scale_fill_distiller(palette = colors, direction = -1, na.value = na_color,
                                      name = "Methylation Level", limits = c(0,1))
    }

  return(g)
}

#' Create a list of GRangesLists for each step in the multiscale plot
#'
#' Output of this function can be passed directly to multiscaleMethylationPlot
#'
#' @param means_dir List of length nSamples. Each element of means.list is itself a list of GRanges of length, nSteps
#' @param samples character vector of samples to include. These samples must match the sample column from samples.tsv in the
#' multiscale plot workflow
#' @param config Config YAML from the multiscale plot workflow
#' @param which GRanges for rtracklayer::import
#'
#' @return a GRangesList with the average values
#'
#' @import GenomicRanges
#'
#' @export
#'
#' @examples
#'
#' # Not a great example, but it works
#' files.loc <- system.file("extdata", package="bisplotti")
#'
#' grp_avg <- multiscaleGroupAverage(
#'     files.loc,
#'     "Heyn_2012_100yr_subset",
#'     paste0(files.loc, "/config_example.yaml"),
#'     NULL
#' )
#'
multiscaleGroupAverage <- function(means_dir, samples, config, which) {
    steps <- .getMultiscaleSteps(config)
    nSteps <- length(steps)

    # Get all of the files for in-group samples from the Snakemake means directory
    files <- unlist(
        lapply(samples,
               function(x) {
                   paste0(means_dir, '/', x, '_', sprintf("%.1f", steps), ".bed.gz")
               }
        )
    )

    # Import each file
    means <- lapply(
        files,
        function(x) {
            return(rtracklayer::import(x, format="bedGraph", which=which))
        }
    )

    # Now split this to a list equal to the length of the samples2grep
    # Each element of means.list is a list of nStep GRanges (similar to if one sample or the data from the stats/
    # directory is imported)
    means.list <- split(means, ceiling(seq_along(means) / nSteps))

    if (!length(samples) == length(means.list)) stop("Length of the means list should be the same as the number of samples.")

    # Apply .getMultiscaleStepGroupAverage for each step size to get a list of GRanges of length nStep
    group.avg.mean <- lapply(seq_along(1:nSteps), .getMultiscaleStepGroupAverage, means.list=means.list)

    stopifnot(length(group.avg.mean) == nSteps)

    group.avg.mean.grl <- as(group.avg.mean, "GRangesList")
    names(group.avg.mean.grl) <- steps

    return(group.avg.mean.grl)
}

# Helpful function: get the bin resolution from the multiscale methylation Snakemake workflow config file
.getMultiscaleSteps <- function(config) {
    if(!grep('\\.yaml$',config)) stop("config must be a yaml file")
    multiscale.config <- yaml::yaml.load_file(config)

    # Get the number of steps
    steps <- seq(multiscale.config$bounds$min, multiscale.config$bounds$max, by=multiscale.config$bounds$step)
    return(steps)
}

# Helper function: get the average for a list of GRanges
.getMultiscaleStepGroupAverage <- function(means.list, i) {
    if (!is(means.list, "list")) stop("means.list needs to be a list")
    if (!is(means.list[[1]], "list")) stop("means.list[[1]] needs to be a list of lists.")
    if (!is(means.list[[1]][[1]], "GRanges")) stop("means.list[[1]][[1]] needs GRanges object.")

    step.gr <- lapply(means.list, `[[`, i) # get the i-th step for all samples

    # Now get the average
    sum <- rep(0,length(step.gr[[1]]))
    for (j in 1:length(step.gr)) {
        sum <- sum + mcols(step.gr[[j]])$score
    }
    avg <- sum / length(step.gr)

    return.gr <- means.list[[1]][[i]] # take the first GRanges for step i
    mcols(return.gr)$score <- avg # replace the score with the average

    return(return.gr)
}
