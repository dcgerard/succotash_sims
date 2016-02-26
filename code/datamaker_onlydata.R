datamaker_counts_only <- function(args) {
    dfargs <- default_datamaker_args(args)

    ## rawdata1 <- readtissue(dfargs$path, dfargs$tissue[1])
    rawdata1 <- read.table(paste0(dfargs$path, dfargs$tissue[1], ".txt"), header = TRUE)   
    if (length(dfargs$tissue) > 1) {
        ## rawdata2 <- readtissue(dfargs$path, dfargs$tissue[2])
        rawdata2 <- read.table(paste0(dfargs$path, dfargs$tissue[2], ".txt"), header = TRUE)

        if (is.null(dfargs$Nsamp)) {
            dfargs$Nsamp <- min(dim(rawdata1)[2], dim(rawdata2)[2])
        }
        if (dim(rawdata1)[2] < dfargs$Nsamp | dim(rawdata2)[2] < dfargs$Nsamp) {
            stop("Not enough samples in the raw dataset!")
        }

        if (dfargs$nullpi == 0) {
            ## All genes are alternatives
            temp1 <- selectsample(rawdata1, dfargs$Nsamp, dfargs$breaksample)
            counts1 <- temp1$counts
            subsample1 <- temp1$subsample
            rm(temp1)
            temp2 <- selectsample(rawdata2, dfargs$Nsamp, dfargs$breaksample)
            counts2 <- temp2$counts
            subsample2 <- temp2$subsample
            rm(temp2)

            counts <- cbind(counts1, counts2)
            subsample <- cbind(subsample1, subsample2)
        } else {
            ## Some genes are nulls, some are alternatives
            temp1 <- selectsample(rawdata1, 2 * dfargs$Nsamp, dfargs$breaksample)
            counts1 <- temp1$counts
            subsample1 <- temp1$subsample
            rm(temp1)
            temp2 <- selectsample(rawdata2, dfargs$Nsamp, dfargs$breaksample)
            counts2 <- temp2$counts
            subsample2 <- temp2$subsample
            rm(temp2)
            counts <- cbind(counts1, counts2)
            subsample <- cbind(subsample1, subsample2)
        }
    } else {
        if (is.null(dfargs$Nsamp)) {
            dfargs$Nsamp <- floor(dim(rawdata1)[2] / 2)
        }
        if (dim(rawdata1)[2] < 2 * dfargs$Nsamp) {
            stop("Not enough samples in the raw dataset!")
        }


        temp <- selectsample(rawdata1, 2 * dfargs$Nsamp, dfargs$breaksample)
        counts <- temp$counts
        subsample <- temp$subsample
        rm(temp)
    }

    ## Remove genes without any reads
    subsample <- subsample[apply(counts, 1, sum) > 0, ]
    counts <- counts[apply(counts, 1, sum) > 0, ]

    ## Take the top Ngene high-expressed genes
    if (!is.null(dfargs$Ngene)) {
        dfargs$Ngene <- min(dfargs$Ngene, dim(counts)[1])
        subsample <- subsample[sort(order(rowSums(counts), decreasing = TRUE)[(dfargs$skip_genes + 1):(dfargs$Ngene + dfargs$skip_genes)]), ]
        counts <- counts[sort(order(rowSums(counts), decreasing = TRUE)[(dfargs$skip_genes + 1):(dfargs$Ngene + dfargs$skip_genes)]), ]
    }
    dfargs$Ngene <- dim(counts)[1]

    ## Model's design: Nsamp samples for group A and Nsamp samples for group B
    condition <- factor(rep(1:2, each = dfargs$Nsamp))
    design <- stats::model.matrix(~condition)

    ## Ground truth of null hypotheses: beta_g=0
    null <- rep(0, dfargs$Ngene)
    null[sample(dfargs$Ngene, round(dfargs$Ngene * dfargs$nullpi))] <- 1
    ## default is that dfargs$nullpi is 1, so all genes are null genes -- dcg

    ## Poisson thinning (optional)
    counts <- pois_thinning(counts, dfargs, null) ## does nothing if args$poisthin == FALSE -- dcg

    ## Mix null and alternative genes from different samples (optional)
    mix <- mix_sample(counts, dfargs, null, subsample) ## only matters if not all genes are null genes -- dcg
    counts <- mix$counts
    subsample <- mix$subsample



    input <- list(counts = counts, condition = condition)
    meta <- list(null = null, dfargs = dfargs)

    return(list(meta = meta,input = input))
}

default_datamaker_args <- function(args) {
    ## poisthin: flag of Poisson thinning
    if (is.null(args$poisthin)) {
        args$poisthin <- FALSE
    }

    ## number of top genes to skip
    if (is.null(args$skip_genes)) {
        args$skip_genes <- 0
    }

    ## log2foldmean, log2foldsd: Poisson thinning params
    if (args$poisthin == TRUE) {
        if (is.null(args$log2foldmean)) {
            args$log2foldmean <- 0
        }
        if (is.null(args$log2foldsd)) {
            args$log2foldsd <- 1
        }
    }

    ## breaksample: flag of each gene randomly select samples
    if (is.null(args$breaksample)) {
        args$breaksample <- FALSE
    }


    ## nullpi: proportion of null genes
    if (is.null(args$nullpi)) {
        if (args$poisthin == TRUE) {
            args$nullpi <- 0.9
        } else if (length(args$tissue) == 1) {
            args$nullpi <- 1
        } else if (length(args$tissue) > 1) {
            args$nullpi <- 0
        } else if (is.null(args$tissue)) {
            args$nullpi <- 1
        }
    }

    ## pseudocounts: add pseudocounts to count matrix
    if (is.null(args$pseudocounts)) {
        args$pseudocounts <- 1
    }

    if(is.null(args$sig_diag)) {
        args$sig_diag <- rep(1, args$Ngene)
    }

    if(is.null(args$sig_alpha)) {
        args$sig_alpha <- 1
    }

    if(is.null(args$beta0)) {
         args$beta0 <- 10
    }

    if(is.null(args$get_null)) {
        args$get_null <- TRUE
    }

    return(args)
}

selectsample <- function(counts, Nsamp, breaksample) {
    if (breaksample == FALSE) {
        subsample <- sample(1:dim(counts)[2], Nsamp)
        counts <- counts[, subsample]
        subsample <- t(matrix(rep(subsample, dim(counts)[1]), ncol = dim(counts)[1]))
    } else {
        temp <- t(apply(counts, 1, sampleingene, Nsamp = Nsamp))
        counts <- temp[, 1:Nsamp]
        subsample <- temp[, (Nsamp + 1):(2 * Nsamp)]
    }
    return(list(counts = counts, subsample = subsample))
}

pois_thinning <- function(counts, args, null) {
    if (args$poisthin == TRUE) {
        log2foldchanges <- rnorm(sum(!null), mean = args$log2foldmean, sd = args$log2foldsd)
        foldchanges <- 2 ^ log2foldchanges

        ## thin group A
        counts[which(!null)[log2foldchanges > 0], 1:args$Nsamp] <-
            matrix(rbinom(sum(log2foldchanges >
            0) * args$Nsamp, size = c(as.matrix(counts[which(!null)[log2foldchanges >
            0], 1:args$Nsamp])), prob = rep(1 / foldchanges[log2foldchanges > 0], args$Nsamp)),
            ncol = args$Nsamp)
        ## thin group B
        counts[which(!null)[log2foldchanges < 0], (args$Nsamp + 1):(2 * args$Nsamp)] <-
            matrix(rbinom(sum(log2foldchanges <
            0) * args$Nsamp, size = c(as.matrix(counts[which(!null)[log2foldchanges <
            0], (args$Nsamp + 1):(2 * args$Nsamp)])), prob = rep(foldchanges[log2foldchanges <
            0], args$Nsamp)), ncol = args$Nsamp)

    }
    return(counts)
}

mix_sample <- function(counts, args, null, subsample) {
    if (args$nullpi < 1 & args$nullpi > 0 & args$breaksample == TRUE) {
        newcounts <- matrix(rep(0, args$Ngene * 2 * args$Nsamp), nrow = args$Ngene)
        newcounts[as.logical(null), ] <- counts[as.logical(null), 1:(2 * args$Nsamp)]
        newcounts[!null, ] <- counts[!null, c(1:args$Nsamp, (2 * args$Nsamp + 1):(3 *
            args$Nsamp))]
        counts <- newcounts
        newsubsample <- matrix(rep(0, args$Ngene * 2 * args$Nsamp), nrow = args$Ngene)
        newsubsample[as.logical(null), ] <- subsample[as.logical(null), 1:(2 * args$Nsamp)]
        newsubsample[!null, ] <- subsample[!null, c(1:args$Nsamp, (2 * args$Nsamp +
            1):(3 * args$Nsamp))]
        subsample <- newsubsample
        rm(newcounts)
        rm(newsubsample)
    }
    return(list(counts = counts, subsample = subsample))
}

fit_ash <- function(out_obj) {
    if(is.null(out_obj$sebetahat)) { ## some methods do not return sebetahat.
        return()
    } else {
        ash_out <- ashr::ash(betahat = out_obj$betahat, sebetahat = out_obj$sebetahat,
                             df = out_obj$df)
        return(ash_out)
    }
}

#' Get wald p-values and adjust by benjamini and hochberg.
fit_freq_methods <- function(out_obj) {
    if (is.null(out_obj$pvalue)) { ## if given pvalues, use those
        if (is.null(out_obj$df)) {
            p_values <- 2 * (1 - pnorm(abs(out_obj$betahat / out_obj$sebetahat)))
        } else {
            p_values <- 2 * (1 - pt(abs(out_obj$betahat / out_obj$sebetahat), df = out_obj$df))
        }
    } else {
        p_values <- out_obj$pvalue
    }
    
    q_bh <- stats::p.adjust(p_values, method = "BH")

    if (all(is.na(p_values))) {
        q_storey <- NULL
    } else {
        q_storey <- qvalue::qvalue(p_values)
    }
    return(list(p_values = p_values, q_bh = q_bh, q_storey = q_storey))
}

extract_ashpi0 <- function(ash_obj) {
    return(ash_obj$fitted.g$pi[1])
}

extract_ashlfdr <- function(ash_obj) {
    return(ash_obj$lfdr)
}

extract_qvaluepi0 <- function(freq_obj) {
    return(freq_obj$q_storey$pi0)
}

extract_pvalues <- function(freq_obj) {
    return(freq_obj$p_values)
}
