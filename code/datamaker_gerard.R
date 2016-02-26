###########################
## Filename: datamaker_gerard.R
## Created by: David Gerard
## Synopsis: Modify Mengyin's code to include many other estimation procedures
###########################

library(edgeR)
library(limma)
library(RUVSeq) 
library(sva)
library(DESeq)
library(DESeq2)
library(succotashr)
library(cate)
library(leapp)






#' Generate (null or alternative) count data from GTEX data and fit a
#' few methods.
#'
#' A function to generate a G*(2N) GTEx count matrix, with G genes and
#' 2N samples for 2 groups.  First N samples for group/condition A,
#' and the rest N samples for group B.
#'
#' @param args is a list of arguments:
#' @param ..path: (Required!) Data path (I put the GTEx data files
#'     under the directory path/gtex/).
#' @param ..tissue: (Required!) String of one tissue name, or a vector
#'     of 2 tissue names.  If inputting one tissue name, then all 2N
#'     samples are selected from this tissue's sample.  If inputting
#'     two tissue names, then first N samples are selected from
#'     tissue1 and last N samples from tissue2 or from a mixture of
#'     tissue1 and tissue2 (when \code{breaksample==TRUE &
#'     0<nullpi<1}).
#' @param ..Nsamp: N, Number of samples in each group (so the total
#'     number of samples is 2*N).  If Nsamp==NULL, select all samples
#'     of the tissue(s).
#' @param ..Ngene: Number of top high-expressed genes that will be
#'     selected. If Ngene==NULL, select all genes.
#' @param ..breaksample: Flag, whether select each gene's counts from
#'     different GTEx samples.  This will break the possible within
#'     sample correlation. The default value is FALSE.
#' @param ..poisthin: Flag, whether use Poisson thinning to change
#'     some counts.  Need specify 2 parameters for thinning:
#'     log2foldmean, log2foldsd, and the proportion of genes for
#'     thinning: nullpi. The default value is FALSE.
#' @param ..log2foldmean: Only works when poisthin==TRUE. Mean
#'     parameter for log2-fold change.
#' @param ..log2foldsd: Only works when poisthin==TRUE. SD parameter
#'     for log2-fold change.
#' @param ..nullpi: Proportion of genes that are nulls.  Only works
#'     when tissue is one tissue name and poisthin==TRUE, or tissue
#'     contains two tissue names and breaksample==TRUE.
#' @param ..pseudocounts: Add pseudocounts to the count matrix. The
#'     default value is 1.
#' @param ..RUV.k: Number of surrogate variables for RUV. The default
#'     value is round(log2(Nsamp)).
#'
#' @return a list of input info and meta info for dscr.
#' 
#' \code{input}: a list of G*2N count matrix (count) and a 2N vector
#'         of group assignment (condition), and a list called
#'         \code{..beta_fits} with at least three elements: estimated
#'         effects (\code{betahat}), sd (\code{sebetahat}) and degree
#'         of freedom (\code{df}) from several different methods.
#'
#' \code{meta}: a list of the subsample index of selected GTEx samples
#'         (subsample), the true null/alternative info for each gene
#'         (null), and the input arguments of datamaker function
#'         (after setting default).
#'
#' @examples Compare 50 Lung samples against 50 Ovary samples (only
#'     select top 10000 highly-expressed genes):
#'     \code{args=list(tissue=c('Lung','Ovary'), Nsamp=50,
#'     Ngene=10000, path='/mnt/lustre/home/mengyin')}
#'
#' More examples in
#' https://github.com/mengyin/dscr-gtex-total/blob/master/scenarios.R.
datamaker <- function(args) {
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

    ## take log of counts
    log_counts <- as.matrix(log(counts + 1))

    ## Voom transformation
    voom_out <- voom_transform(counts, condition)

    ## Quasi-binomial glm
    qb_out <- quasi_binom(counts, condition)

    ## Myrna (Langmead et al. '10) & Quasi-binomial glm Use log(75th quantile of
    ## samples' counts) as covariate
    W.Myrna <- apply(counts, 2, function(x) log(quantile(x[x > 0], 0.75)))
    Myrnaqb_out <- quasi_binom(counts, condition, W = W.Myrna)

    ## Use log(75th quantile of samples' counts) as offsets
    offset.Myrnaoff <- apply(counts, 2, function(x) quantile(x[x > 0], 0.75))
    Myrnaoffqb_out <- quasi_binom(counts, condition, W = NULL, offset = offset.Myrnaoff)

    ## estimate number of confounders using log-counts
    dfargs$num.sv <- sva::num.sv(log_counts, mod = stats::model.matrix(~condition), method = "be")

    ## RUV & voom
    halfnull <- rep(0, length(null))  ## Use half of the true nulls to do supervised RUV/SVA
    halfnull[which(null == 1)[1:floor(length(which(null == 1)) / 2)]] <- 1
    W.RUV <- RUV_factor(counts, dfargs, halfnull) 
    RUVvoom_out <- voom_transform(counts, condition, W = W.RUV)

    ## SVA & voom
    W.SVA <- SVA_factor(counts, condition, dfargs)
    SVAvoom_out <- voom_transform(counts, condition, W = W.SVA)

    ## RUV & quasi-binomial glm
    RUVqb_out <- quasi_binom(counts, condition, W = W.RUV)

    ## SVA & quasi-binomial glm
    SVAqb_out <- quasi_binom(counts, condition, W = W.SVA)

    ## Get sebetahat from edgeR.glm (infer from betahat & pval)
    #edgeRglm_out <- edgeR_glmest(counts, condition, dfargs)

    ## Get sebetahat from DESeq2 (infer from betahat & pval)
    DESeq2glm_out <- DESeq2_glmest(counts, condition, dfargs)

    ## meta data
    meta <- list(null = null, dfargs = dfargs)

    ## ols methods
    ols_out <- get_ols(log_counts, condition)

    design.ruv <- cbind(as.numeric(condition) - 1, W.RUV)
    ols_ruv_out <- get_ols(log_counts, design.ruv)

    design.sva <- cbind(as.numeric(condition) - 1, W.SVA)
    ols_sva_out <- get_ols(log_counts, design.sva)

    ## robust regression version of CATE on log(counts + 1)
    cate_rr <- cate::cate(~trt, Y = t(log_counts),
                          X.data = data.frame(trt = as.numeric(condition) - 1),
                          r = dfargs$num.sv, fa.method = "ml", adj.method = "rr")
    cate_rr_out <- list()
    cate_rr_out$betahat <- cate_rr$beta
    cate_rr_out$sebetahat <- sqrt(cate_rr$beta.cov.row * cate_rr$beta.cov.col) / sqrt(dfargs$Nsamp)
    cate_rr_out$pvalue <- cate_rr$beta.p.value ## median adjusted to correct for standard deviation
                                               ## being biased small

    ## Control gene version of CATE on log(counts + 1)
    cate_nc <- cate::cate(~trt, Y = t(log_counts),
                          X.data = data.frame(trt = as.numeric(condition) - 1),
                          r = dfargs$num.sv, fa.method = "ml", adj.method = "nc", nc = as.logical(halfnull))
    cate_nc_out <- list()
    cate_nc_out$betahat <- cate_nc$beta
    cate_nc_out$sebetahat <- sqrt(cate_nc$beta.cov.row * cate_nc$beta.cov.col) / sqrt(dfargs$Nsamp)
    cate_nc_out$pvalue <- cate_nc$beta.p.value ## median adjusted to correct for standard deviation
                                               ## being biased small

    ## LEAPP this takes the longest out of all of the methods in here, especially for small p
    ## (probably convergence issues)
    leapp_sparse <- leapp::leapp(dat = log_counts, pred.prim = as.numeric(condition) - 1,
                                 num.fac = dfargs$num.sv, method = "soft")
    leapp_sparse_out <- list()
    leapp_sparse_out$betahat <- leapp_sparse$gamma
    leapp_sparse_out$pvalue <- leapp_sparse$p
    
    leapp_ridge <- leapp::leapp(dat = log_counts, pred.prim = as.numeric(condition) - 1,
                                num.fac = dfargs$num.sv, method = "hard", sparse = FALSE)
    leapp_ridge_out <- list()
    leapp_ridge_out$betahat <- leapp_ridge$gamma
    leapp_ridge_out$pvalue <- leapp_ridge$p
    
    
    ## initialize list of fits
    beta_fits <- list(voom_out = voom_out, qb_out = qb_out, Myrnaqb_out = Myrnaqb_out,
                      Myrnaoffqb_out = Myrnaoffqb_out, RUVvoom_out = RUVvoom_out,
                      SVAvoom_out = SVAvoom_out, RUVqb_out = RUVqb_out, SVAqb_out = SVAqb_out,
                      DESeq2glm_out = DESeq2glm_out,
                      ols_out = ols_out, ols_ruv_out = ols_ruv_out, ols_sva_out = ols_sva_out,
                      cate_rr_out = cate_rr_out, cate_nc_out = cate_nc_out,
                      leapp_sparse_out = leapp_sparse_out, leapp_ridge_out = leapp_ridge_out)

    ## input data
    input <- list(counts = counts, condition = condition, v = voom_out$v, W.RUV = W.RUV,
                  W.SVA = W.SVA, W.Myrna = W.Myrna)

    data <- list(meta = meta, input = input, beta_fits = beta_fits)
    return(data)
}

## Set default arguments for datamaker function
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

## Poisson thinning
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

## Mix null and alternative genes from different samples
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

## Voom transformation
voom_transform <- function(counts, condition, W = NULL) {
    dgecounts <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts, group = condition))

    if (is.null(W)) {
        design <- model.matrix(~condition)
    } else {
        design <- model.matrix(~condition + W)
    }

    v <- limma::voom(dgecounts, design, plot = FALSE)
    lim <- limma::lmFit(v)
    ## zdat.voom <- apply(cbind(v$E,v$weights),1,wls.wrapper,g=condition) betahat.voom
    ## <- zdat.voom[1,] sebetahat.voom = zdat.voom[2,]
    betahat.voom <- lim$coefficients[, 2]
    sebetahat.voom <- lim$stdev.unscaled[, 2] * lim$sigma
    df.voom <- length(condition) - 2 - (!is.null(W))

    return(list(betahat = betahat.voom, sebetahat = sebetahat.voom, df = df.voom,
        v = v))
}

## Weighted least squares regression g: formula ynweights: matrix of response y
## and corresponding weights
wls.wrapper <- function(ynweights, g, ...) {
    y <- ynweights[1:(length(ynweights) / 2)]
    weights <- ynweights[(length(ynweights) / 2 + 1):length(ynweights)]
    y.wls <- lm(y ~ g, weights = weights, ...)

    ## slope estimate & standard error
    c <- as.vector(t(summary(y.wls)$coeff[2, 1:2]))
    return(c)
}

## Quasi-binomial glm
quasi_binom <- function(counts, condition, W = NULL, offset = NULL) {
    zdat.qb <- counts.associate(counts, condition, W = W, offset = offset)
    betahat <- zdat.qb[3, ]
    sebetahat <- zdat.qb[4, ]
    dispersion <- zdat.qb[5, ]
    df <- length(condition) - 2 - (!is.null(W))
    return(list(betahat = betahat, sebetahat = sebetahat, df = df, dispersion = dispersion))
}

## counts is a ngene (or nwindow) by nsample matrix of counts (eg RNAseq) g is an
## nsample vector of group memberships/covariate looks for association between
## rows and covariate
counts.associate <- function(counts, g, W = NULL, offset = NULL, pseudocounts = 1) {
    y.counts <- t(as.matrix(counts))
    col.sum <- apply(y.counts, 2, sum)
    y.counts <- y.counts[, col.sum > 0]  #remove 0 columns
    y.counts <- y.counts + pseudocounts
    if (is.null(offset)) {
        offset <- apply(y.counts, 1, sum)
    }

    y.prop <- y.counts / apply(y.counts, 1, sum)  ## compute proportions
    zdat <- apply(y.prop, 2, glm.binomial.wrapper, g = g, W = W, weights = offset,
        epsilon = 1e-06)  #estimate effect sizes and standard errors
    ## zdat.ash <- ash(zdat[3,],zdat[4,],df=2,method='fdr') #shrink the estimated
    ## effects return(list(zdat=zdat,zdat.ash=zdat.ash))
    return(zdat)
}

glm.binomial.wrapper <- function(y, g, W = NULL, ...) {
    if (is.null(W)) {
        y.glm <- safe.quasibinomial.glm(y ~ g, ...)
    } else {
        y.glm <- safe.quasibinomial.glm(y ~ g + W, ...)
    }
    return(c(get.coeff(y.glm), summary(y.glm)$dispersion))
}


## fill NAs with 0s (or other value)
fill.nas <- function(x, t = 0) {
    x[is.na(x)] <- t
    return(x)
}

## get estimates and standard errors from a glm object return NAs if not converged
get.coeff <- function(x.glm) {
    c <- as.vector(t(summary(x.glm)$coeff[, 1:2]))
    if (x.glm$conv) {
        return(c)
    } else {
        return(rep(NA, length(c)))
    }
}

## use glm to fit quasibinomial, but don't allow for underdispersion!
safe.quasibinomial.glm <- function(formula, forcebin = FALSE, ...) {
    if (forcebin) {
        fm <- glm(formula, family = binomial, ...)
    } else {
        fm <- glm(formula, family = quasibinomial, ...)
        if (is.na(summary(fm)$dispersion) | summary(fm)$dispersion < 1) {
            fm <- glm(formula, family = binomial, ...)
        }
    }
    return(fm)
}

## randomly subsample data for each gene gene: a vector of reads for one gene
## Nsamp: ## of samples wanted
sampleingene <- function(gene, Nsamp) {
    sample <- sample(length(gene), Nsamp)
    return(c(gene[sample], sample))
}

## Randomly select samples counts: full count matrix Nsamp: ## of samples wanted
## breaksample: flag, if select different samples for each gene
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

## Use RUV to estimate confounding factor
## I changed this to use the
## ruv package instead of RUVSeq because RUVSeq requires loading other
## libraries to use individual functions, which is dumb. RUVSeq is
## just RUV2 using the log-counts, so whatever. --- dcg
## RUV_factor <- function(counts_current, condition, args, null) {
##     if (sum(null) > 0) {
##         ruv_out <- ruv::RUV2(Y = t(log(counts_current)), X = matrix(as.numeric(condition) - 1),
##                              ctl = as.logical(null), k = args$num.sv)
        
##         return(W = ruv_out$W)
##     } else {
##         return(W = NULL)
##     }
## }

## Use RUV to estimate confounding factor
RUV_factor = function(counts, args, null) {
    seq = EDASeq::newSeqExpressionSet(as.matrix(counts[as.logical(null), ]))
    if (sum(null) > 0) {
        controls = rownames(seq)
        differences = matrix(data = c(1:args$Nsamp, (args$Nsamp + 1):(2 * args$Nsamp)), 
            byrow = TRUE, nrow = 2)
        seqRUV = RUVSeq::RUVg(seq, controls, k = args$num.sv)
        return(W = as.matrix(Biobase::pData(seqRUV)))
    } else {
        return(W = NULL)
    }
}

## Use SVA to estimate confounding factor
## changed nullg so default is not to have control genes
SVA_factor <- function(counts, condition, args, nullg = NULL) {
    mod1 <- model.matrix(~condition)
    mod0 <- cbind(mod1[, 1])
    
    if (!is.null(nullg)) {
        svseq_out <- sva::svaseq(counts, mod1, mod0, control = nullg, n.sv = args$num.sv)
    } else {
        svseq_out <- sva::svaseq(as.matrix(counts), mod1, mod0, n.sv = args$num.sv)
    }

    if (svseq_out$n.sv > 0) {
        return(W = svseq_out$sv)
    } else {
        return(W = NULL)
    }
}

## Get sebetahat from edgeR.glm (infer from betahat & pval)
edgeR_glmest <- function(counts, condition, args) {
    design <- model.matrix(~condition)
    y <- edgeR::DGEList(counts = counts + args$pseudocounts, group = condition)
    y <- edgeR::calcNormFactors(y)
    y <- edgeR::estimateGLMCommonDisp(y, design)
    y <- edgeR::estimateGLMTrendedDisp(y, design)
    y <- edgeR::estimateGLMTagwiseDisp(y, design)
    fit <- edgeR::glmFit(y, design)
    lrt <- edgeR::glmLRT(fit, coef = 2)
    betahat.edgeRglm <- fit$coef[, 2]
    df.edgeRglm <- length(condition) - 2
    tscore <- qt(1 - lrt$table$PValue / 2, df = df.edgeRglm)
    sebetahat.edgeRglm <- abs(betahat.edgeRglm / tscore)
    return(list(betahat = betahat.edgeRglm, sebetahat = sebetahat.edgeRglm, df = df.edgeRglm))
}

## Get sebetahat from DESeq2 (infer from betahat & pval)
DESeq2_glmest <- function(counts, condition, args) {
    cond <- condition
    dds <- DESeq2::DESeqDataSetFromMatrix(counts + args$pseudocounts,
                                         S4Vectors::DataFrame(cond), ~cond)
    dds <- DESeq2::estimateSizeFactors(dds)
    dds <- DESeq2::estimateDispersions(dds, fitType = "local")
    dds <- DESeq2::nbinomWaldTest(dds)

    #mcols(mcols(dds, use.names = TRUE))
    #coef(dds)
    #dout <- DESeq(dds)
    #dcoef <- coef(dout, SE = TRUE)
    #class(dcoef)
    res <- DESeq2::results(dds, cooksCutoff = FALSE)
    betahat.DESeq2 <- res$log2FoldChange
    sebetahat.DESeq2 <- res$lfcSE
    return(list(betahat = betahat.DESeq2, sebetahat = sebetahat.DESeq2))
}


## Extract dataset for a specific tissue from the GTEx reads txt file Note: use
## GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct_new.txt, which removes
## the first 3 lines from
## GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.txt
readtissue <- function(path, tissue) {
    tis <- read.table(paste0(path, "/data/sample_tissue.txt"), header = TRUE)

    tissue.idx <- grep(tissue, tis[, 2], fixed = TRUE)
    cols <- rep("NULL", dim(tis)[1])
    cols[tissue.idx] <- "numeric"

    data <- read.table(paste0(path, "/gtex/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct_new.txt"), 
        colClasses <- cols, header = TRUE)
    return(data)
}

##' Simple ordinary least squares.
##'
##' @param log_counts A matrix of numerics. The responses.
##' @param condition A matrix of numerics. The predictors.
get_ols <- function(log_counts, condition) {
    limma_out <- limma::lmFit(log_counts, model.matrix(~condition))
    betahat <- limma_out$coefficients[, 2]
    sebetahat <- limma_out$stdev.unscaled[, 2] * limma_out$sigma
    df <- limma_out$df.residual
    return(list(betahat = betahat, sebetahat = sebetahat, df = df))
}

#' Fit a bunch of ash procedures.
#'
#' @param out_obj A list with at least three objects.
#' @param ..betahat A vector of numerics. The fits of beta.
#' @param ..sebetahat A vector of numerics. The fits of the standard errors of \code{betahat}s.
#' @param ..df A numeric or vector of numerics. The degrees of freedom of \code{betahat}.
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
    } else if(max(p_values, na.rm = TRUE) < 0.9) {
        ## need to do this so that don't get thrown an error
        q_storey <- qvalue::qvalue(p_values, lambda = quantile(p_values, 0.75))
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


### test data, remove when done
## args <- list()
## args$tissue <- "Liver"
## args$Nsamp <- 20
## args$path <- "../data/"
## args$Ngene <- 100


## data_out <- datamaker(args)
## ash_out <- lapply(data_out$beta_fits, fit_ash)
## sapply(ash_out, extract_ashpi0)
## sapply(ash_out, extract_ashlfdr)
## freq_out <- lapply(data_out$beta_fits, fit_freq_methods)
## sapply(freq_out, extract_qvaluepi0)

## succ_out <- succotashr::succotash(Y = t(as.matrix(log(data_out$input$counts))),
##                      X = model.matrix(~data_out$input$condition),
##                       k = data_out$meta$dfargs$num.sv,
##                       fa_method = "quasi_mle", lambda0 = 10)
## plot(succ_out$tau_seq, succ_out$pi, type = "h")
