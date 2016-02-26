

library(edgeR)
library(limma)
library(RUVSeq)
library(sva)
library(DESeq)
library(DESeq2)
library(succotashr)



#' Generate (null or alternative) count data from GTEX data and fit a few methods.
#'
#' A function to generate a G*(2N) GTEx count matrix, with G genes and
#' 2N samples for 2 groups.  First N samples for group/condition A, and the rest
#' N samples for group B.
#'
#' @param args is a list of arguments:
#' @param ..path: (Required!) Data path (I put the GTEx data files under the directory 
#'         path/gtex/).
#' @param ..tissue: (Required!) String of one tissue name, or a vector of 2 tissue names.
#'         If inputting one tissue name, then all 2N samples are selected from this 
#'         tissue's sample.  If inputting two tissue names, then first N samples are 
#'         selected from tissue1 and last N samples from tissue2 or from a mixture of 
#'         tissue1 and tissue2 (when \code{breaksample==TRUE & 0<nullpi<1}).
#' @param ..Nsamp: N, Number of samples in each group (so the total number of samples is 
#'         2*N).  If Nsamp==NULL, select all samples of the tissue(s).
#' @param ..Ngene: Number of top high-expressed genes that will be selected. If 
#'         Ngene==NULL, select all genes.
#' @param ..breaksample: Flag, whether select each gene's counts from different GTEx 
#'         samples.  This will break the possible within sample correlation. The default 
#'         value is FALSE.
#' @param ..poisthin: Flag, whether use Poisson thinning to change some counts.  Need 
#'         specify 2 parameters for thinning: log2foldmean, log2foldsd, and the proportion
#'         of genes for thinning: nullpi. The default value is FALSE.
#' @param ..log2foldmean: Only works when poisthin==TRUE. Mean parameter for log2-fold 
#'         change.
#' @param ..log2foldsd: Only works when poisthin==TRUE. SD parameter for log2-fold 
#'         change.
#' @param ..nullpi: Proportion of genes that are nulls.  Only works when tissue is one 
#'         tissue name and poisthin==TRUE, or tissue contains two tissue names and 
#'         breaksample==TRUE.
#' @param ..pseudocounts: Add pseudocounts to the count matrix. The default value is 1.
#' @param ..RUV.k: Number of surrogate variables for RUV. The default value is 
#'         round(log2(Nsamp)).
#'
#' @return a list of input info and meta info for dscr.
#' 
#' \code{input}: a list of G*2N count matrix (count) and a 2N vector of group 
#'         assignment (condition), and estimated effects (betahat.XXX), sd (sebetahat.XXX)
#'         and degree of freedom (df.XXX) from several different methods.
#'
#' \code{meta}: a list of the subsample index of selected GTEx samples (subsample), the
#'         true null/alternative info for each gene (null), and the input arguments of 
#'         datamaker function (after setting default).
#'
#' @examples
#' Compare 50 Lung samples against 50 Ovary samples (only select top 
#'         10000 highly-expressed genes): \code{args=list(tissue=c('Lung','Ovary'), Nsamp=50, Ngene=10000, path='/mnt/lustre/home/mengyin')}
#'
#' More examples in https://github.com/mengyin/dscr-gtex-total/blob/master/scenarios.R.
datamaker <- function(args) {
    dfargs = default_datamaker_args(args)

    ## rawdata1 = readtissue(dfargs$path, dfargs$tissue[1])
    rawdata1 = read.table(paste0(dfargs$path, dfargs$tissue[1], ".txt"), header = TRUE)
    
    if (length(dfargs$tissue) > 1) {
        ## rawdata2 = readtissue(dfargs$path, dfargs$tissue[2])
        rawdata2 = read.table(paste0(dfargs$path, dfargs$tissue[2], ".txt"), header = TRUE)
        
        if (is.null(dfargs$Nsamp)) {
            dfargs$Nsamp = min(dim(rawdata1)[2], dim(rawdata2)[2])
        }
        if (dim(rawdata1)[2] < dfargs$Nsamp | dim(rawdata2)[2] < dfargs$Nsamp) {
            stop("Not enough samples in the raw dataset!")
        }
        
        
        if (dfargs$nullpi == 0) {
            ## All genes are alternatives
            temp1 = selectsample(rawdata1, dfargs$Nsamp, dfargs$breaksample)
            counts1 = temp1$counts
            subsample1 = temp1$subsample
            rm(temp1)
            temp2 = selectsample(rawdata2, dfargs$Nsamp, dfargs$breaksample)
            counts2 = temp2$counts
            subsample2 = temp2$subsample
            rm(temp2)
            
            counts = cbind(counts1, counts2)
            subsample = cbind(subsample1, subsample2)
        } else {
            ## Some genes are nulls, some are alternatives
            temp1 = selectsample(rawdata1, 2 * dfargs$Nsamp, dfargs$breaksample)
            counts1 = temp1$counts
            subsample1 = temp1$subsample
            rm(temp1)
            temp2 = selectsample(rawdata2, dfargs$Nsamp, dfargs$breaksample)
            counts2 = temp2$counts
            subsample2 = temp2$subsample
            rm(temp2)
            counts = cbind(counts1, counts2)
            subsample = cbind(subsample1, subsample2)
        }
    } else {
        if (is.null(dfargs$Nsamp)) {
            dfargs$Nsamp = floor(dim(rawdata1)[2]/2)
        }
        if (dim(rawdata1)[2] < 2 * dfargs$Nsamp) {
            stop("Not enough samples in the raw dataset!")
        }
        
        
        temp = selectsample(rawdata1, 2 * dfargs$Nsamp, dfargs$breaksample)
        counts = temp$counts
        subsample = temp$subsample
        rm(temp)
    }
    
    ## Remove genes without any reads
    subsample = subsample[apply(counts, 1, sum) > 0, ]
    counts = counts[apply(counts, 1, sum) > 0, ]
    
    ## Take the top Ngene high-expressed genes
    if (!is.null(dfargs$Ngene)) {
        dfargs$Ngene = min(dfargs$Ngene, dim(counts)[1])
        subsample = subsample[sort(order(rowSums(counts), decreasing = TRUE)[1:dfargs$Ngene]), 
                              ]
        counts = counts[sort(order(rowSums(counts), decreasing = TRUE)[1:dfargs$Ngene]), 
                        ]
    }
    dfargs$Ngene = dim(counts)[1]
    
    ## Model's design: Nsamp samples for group A and Nsamp samples for group B
    condition = factor(rep(1:2, each = dfargs$Nsamp))
    design = model.matrix(~condition)
    
    ## Ground truth of null hypotheses: beta_g=0
    null = rep(0, dfargs$Ngene)
    null[sample(dfargs$Ngene, round(dfargs$Ngene * dfargs$nullpi))] = 1 ## default is that dfargs$nullpi is 1, so all genes are null genes -- dcg 
    
    ## Poisson thinning (optional)
    counts = pois_thinning(counts, dfargs, null) ## does nothing if args$poisthin == FALSE -- dcg
    
    ## Mix null and alternative genes from different samples (optional)
    mix = mix_sample(counts, dfargs, null, subsample) ## only matters if not all genes are null genes -- dcg
    counts = mix$counts
    subsample = mix$subsample
    
    ## Voom transformation
    voom = voom_transform(counts, condition)
    
    # Quasi-binomial glm
    qb = quasi_binom(counts, condition)
    
    # Myrna (Langmead et al. '10) & Quasi-binomial glm Use log(75th quantile of
    # samples' counts) as covariate
    W.Myrna = apply(counts, 2, function(x) log(quantile(x[x > 0], 0.75)))
    Myrnaqb = quasi_binom(counts, condition, W = W.Myrna)
    # Use log(75th quantile of samples' counts) as offsets
    offset.Myrnaoff = apply(counts, 2, function(x) quantile(x[x > 0], 0.75))
    Myrnaoffqb = quasi_binom(counts, condition, W = NULL, offset = offset.Myrnaoff)
    
    # RUV & voom
    halfnull = rep(0, length(null))  # Use half of the true nulls to do supervised RUV/SVA
    halfnull[which(null == 1)[1:floor(length(which(null == 1))/2)]] = 1
    W.RUV = RUV_factor(counts, dfargs, halfnull) ## only uses first confounder
    RUVvoom = voom_transform(counts, condition, W = W.RUV)
    
    # SVA & voom
    W.SVA = SVA_factor(counts, condition, dfargs, halfnull)
    SVAvoom = voom_transform(counts, condition, W = W.SVA)
    
    # RUV & quasi-binomial glm
    RUVqb = quasi_binom(counts, condition, W = W.RUV)
    
    # SVA & quasi-binomial glm
    SVAqb = quasi_binom(counts, condition, W = W.SVA)
    
    # Get sebetahat from edgeR.glm (infer from betahat & pval)
    edgeRglm = edgeR_glmest(counts, condition, dfargs)
    
    # Get sebetahat from DESeq2 (infer from betahat & pval)
    DESeq2glm = DESeq2_glmest(counts, condition, dfargs)
    
    # meta data
    meta = list(null = null, dfargs = dfargs)
    
    # input data
    input = list(counts = counts, condition = condition, v = voom$v, betahat.voom = voom$betahat, 
        sebetahat.voom = voom$sebetahat, df.voom = voom$df, betahat.RUVvoom = RUVvoom$betahat, 
        sebetahat.RUVvoom = RUVvoom$sebetahat, df.RUVvoom = RUVvoom$df, W.RUV = W.RUV, 
        betahat.SVAvoom = SVAvoom$betahat, sebetahat.SVAvoom = SVAvoom$sebetahat, 
        df.SVAvoom = SVAvoom$df, W.SVA = W.SVA, betahat.qb = qb$betahat, sebetahat.qb = qb$sebetahat, 
        df.qb = qb$df, dispersion.qb = qb$dispersion, betahat.RUVqb = RUVqb$betahat, 
        sebetahat.RUVqb = RUVqb$sebetahat, dispersion.RUVqb = RUVqb$dispersion, df.RUVqb = RUVqb$df, 
        W.RUV = W.RUV, betahat.SVAqb = SVAqb$betahat, sebetahat.SVAqb = SVAqb$sebetahat, 
        dispersion.SVAqb = SVAqb$dispersion, df.SVAqb = SVAqb$df, W.SVA = W.SVA, 
        betahat.Myrnaqb = Myrnaqb$betahat, sebetahat.Myrnaqb = Myrnaqb$sebetahat, 
        dispersion.Myrnaqb = Myrnaqb$dispersion, df.Myrnaqb = Myrnaqb$df, W.Myrna = W.Myrna, 
        betahat.Myrnaoffqb = Myrnaoffqb$betahat, sebetahat.Myrnaoffqb = Myrnaoffqb$sebetahat, 
        dispersion.Myrnaoffqb = Myrnaoffqb$dispersion, df.Myrnaoffqb = Myrnaoffqb$df, 
        offset.Myrnaoff = offset.Myrnaoff, betahat.edgeRglm = edgeRglm$betahat, sebetahat.edgeRglm = edgeRglm$sebetahat, 
        df.edgeRglm = edgeRglm$df, betahat.DESeq2glm = DESeq2glm$betahat, sebetahat.DESeq2glm = DESeq2glm$sebetahat, 
        df.DESeq2glm = DESeq2glm$df)
    
    data = list(meta = meta, input = input)
    return(data)
}

#' Same as \code{datamaker} but generate data \code{counts} such that \code{log(counts + 1)} is normal.
#'
#' @param args a list of arguments.
#' @param Nsamp
#' @param Ngene
#' @param Nconfounders
#' @param sig_alpha A positive numeric. The variances of the coefficients of the confounders.
#'                  Default is 1
#' @param sig_diag A vector of positive numerics of length Ngene. The variances of the columns
#'                 (genes). Default is a vector of ones.
#' @param beta0 A positive numeric. The intercept for all regressions.
#' @param get_null A logical. Should we simulate under all null (TRUE) or some alternative (FALSE)
#' @param tau_seq Vector of positive numerics. The mixture standard devations. First element should be 0.
#' @param pi_seq Vector of positive numerics that sum to 1. The mixture proportions.
theo_datamaker <- function(args)
{
    dfargs = default_datamaker_args(args)

    theo_out <- theo_count_dat(dfargs$Nsamp, dfargs$Ngene, dfargs$Nconfounders, dfargs$sig_diag,
                   dfargs$sig_alpha, beta0 = dfargs$beta0, get_null = dfargs$get_null,
                   tau_seq = dfargs$tau_seq, pi_seq = dfargs$pi_seq)

    counts <- t(theo_out$counts)
    condition = factor(theo_out$X[,2])
    design = theo_out$X

    null = (theo_out$beta[2,] < 10^-3) * 1 ## 1 if null gene, 0 otherwise
    

    ## Poisson thinning (optional)
    counts = pois_thinning(counts, dfargs, null) ## does nothing if args$poisthin == FALSE -- dcg
    
    ## Voom transformation
    voom = voom_transform(counts, condition)
    
    # RUV & voom
    halfnull = rep(0, length(null))  # Use half of the true nulls to do supervised RUV/SVA
    halfnull[which(null == 1)[1:floor(length(which(null == 1))/2)]] = 1
    W.RUV = RUV_factor(counts, dfargs, halfnull) ## only uses first confounder
    RUVvoom = voom_transform(counts, condition, W = W.RUV)
    
    # SVA & voom
    W.SVA = SVA_factor(counts, condition, dfargs, halfnull)
    SVAvoom = voom_transform(counts, condition, W = W.SVA)
    
    # Get sebetahat from edgeR.glm (infer from betahat & pval)
    edgeRglm = edgeR_glmest(counts, condition, dfargs)
    
    ## meta data
    meta = list(null = null, dfargs = dfargs, theo_out = theo_out)
    
    # input data
    input = list(counts = counts, condition = condition, v = voom$v, betahat.voom = voom$betahat, 
        sebetahat.voom = voom$sebetahat, df.voom = voom$df, betahat.RUVvoom = RUVvoom$betahat, 
        sebetahat.RUVvoom = RUVvoom$sebetahat, df.RUVvoom = RUVvoom$df, W.RUV = W.RUV, 
        betahat.SVAvoom = SVAvoom$betahat, sebetahat.SVAvoom = SVAvoom$sebetahat, 
        df.SVAvoom = SVAvoom$df, W.SVA = W.SVA, 
        W.RUV = W.RUV, W.SVA = W.SVA, betahat.edgeRglm = edgeRglm$betahat,
        sebetahat.edgeRglm = edgeRglm$sebetahat, 
        df.edgeRglm = edgeRglm$df)
    
    data = list(meta = meta, input = input)
    return(data)
}

# Set default arguments for datamaker function
default_datamaker_args = function(args) {
    # poisthin: flag of Poisson thinning
    if (is.null(args$poisthin)) {
        args$poisthin = FALSE
    }
    
    # log2foldmean, log2foldsd: Poisson thinning params
    if (args$poisthin == TRUE) {
        if (is.null(args$log2foldmean)) {
            args$log2foldmean = 0
        }
        if (is.null(args$log2foldsd)) {
            args$log2foldsd = 1
        }
    }
    
    # breaksample: flag of each gene randomly select samples
    if (is.null(args$breaksample)) {
        args$breaksample = FALSE
    }
    
    
    # nullpi: proportion of null genes
    if (is.null(args$nullpi)) {
        if (args$poisthin == TRUE) {
            args$nullpi = 0.9
        } else if (length(args$tissue) == 1) {
            args$nullpi = 1
        } else if (length(args$tissue) > 1) {
            args$nullpi = 0
        } else if (is.null(args$tissue)) {
            args$nullpi = 1
        }
    }
    
    # RUV.k: number of surrogate variables for RUV.
    if (is.null(args$RUV.k)) {
        args$RUV.k = round(log2(args$Nsamp))
    }
    
    # pseudocounts: add pseudocounts to count matrix
    if (is.null(args$pseudocounts)) {
        args$pseudocounts = 1
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

# Poisson thinning
pois_thinning = function(counts, args, null) {
    if (args$poisthin == TRUE) {
        log2foldchanges = rnorm(sum(!null), mean = args$log2foldmean, sd = args$log2foldsd)
        foldchanges = 2^log2foldchanges
        
        # thin group A
        counts[which(!null)[log2foldchanges > 0], 1:args$Nsamp] = matrix(rbinom(sum(log2foldchanges > 
            0) * args$Nsamp, size = c(as.matrix(counts[which(!null)[log2foldchanges > 
            0], 1:args$Nsamp])), prob = rep(1/foldchanges[log2foldchanges > 0], args$Nsamp)), 
            ncol = args$Nsamp)
        # thin group B
        counts[which(!null)[log2foldchanges < 0], (args$Nsamp + 1):(2 * args$Nsamp)] = matrix(rbinom(sum(log2foldchanges < 
            0) * args$Nsamp, size = c(as.matrix(counts[which(!null)[log2foldchanges < 
            0], (args$Nsamp + 1):(2 * args$Nsamp)])), prob = rep(foldchanges[log2foldchanges < 
            0], args$Nsamp)), ncol = args$Nsamp)
        
    }
    return(counts)
}

# Mix null and alternative genes from different samples
mix_sample = function(counts, args, null, subsample) {
    if (args$nullpi < 1 & args$nullpi > 0 & args$breaksample == TRUE) {
        newcounts = matrix(rep(0, args$Ngene * 2 * args$Nsamp), nrow = args$Ngene)
        newcounts[as.logical(null), ] = counts[as.logical(null), 1:(2 * args$Nsamp)]
        newcounts[!null, ] = counts[!null, c(1:args$Nsamp, (2 * args$Nsamp + 1):(3 * 
            args$Nsamp))]
        counts = newcounts
        newsubsample = matrix(rep(0, args$Ngene * 2 * args$Nsamp), nrow = args$Ngene)
        newsubsample[as.logical(null), ] = subsample[as.logical(null), 1:(2 * args$Nsamp)]
        newsubsample[!null, ] = subsample[!null, c(1:args$Nsamp, (2 * args$Nsamp + 
            1):(3 * args$Nsamp))]
        subsample = newsubsample
        rm(newcounts)
        rm(newsubsample)
    }
    return(list(counts = counts, subsample = subsample))
}

# Voom transformation
voom_transform = function(counts, condition, W = NULL) {
    dgecounts = calcNormFactors(DGEList(counts = counts, group = condition))
    
    if (is.null(W)) {
        design = model.matrix(~condition)
    } else {
        design = model.matrix(~condition + W)
    }
    
    v = voom(dgecounts, design, plot = FALSE)
    lim = lmFit(v)
    # zdat.voom = apply(cbind(v$E,v$weights),1,wls.wrapper,g=condition) betahat.voom
    # = zdat.voom[1,] sebetahat.voom = zdat.voom[2,]
    betahat.voom = lim$coefficients[, 2]
    sebetahat.voom = lim$stdev.unscaled[, 2] * lim$sigma
    df.voom = length(condition) - 2 - (!is.null(W))
    
    return(list(betahat = betahat.voom, sebetahat = sebetahat.voom, df = df.voom, 
        v = v))
}

# Weighted least squares regression g: formula ynweights: matrix of response y
# and corresponding weights
wls.wrapper = function(ynweights, g, ...) {
    y = ynweights[1:(length(ynweights)/2)]
    weights = ynweights[(length(ynweights)/2 + 1):length(ynweights)]
    y.wls = lm(y ~ g, weights = weights, ...)
    
    # slope estimate & standard error
    c = as.vector(t(summary(y.wls)$coeff[2, 1:2]))
    return(c)
}

# Quasi-binomial glm
quasi_binom = function(counts, condition, W = NULL, offset = NULL) {
    zdat.qb = counts.associate(counts, condition, W = W, offset = offset)
    betahat = zdat.qb[3, ]
    sebetahat = zdat.qb[4, ]
    dispersion = zdat.qb[5, ]
    df = length(condition) - 2 - (!is.null(W))
    return(list(betahat = betahat, sebetahat = sebetahat, df = df, dispersion = dispersion))
}

# counts is a ngene (or nwindow) by nsample matrix of counts (eg RNAseq) g is an
# nsample vector of group memberships/covariate looks for association between
# rows and covariate
counts.associate = function(counts, g, W = NULL, offset = NULL, pseudocounts = 1) {
    y.counts = t(as.matrix(counts))
    col.sum = apply(y.counts, 2, sum)
    y.counts = y.counts[, col.sum > 0]  #remove 0 columns
    y.counts = y.counts + pseudocounts
    if (is.null(offset)) {
        offset = apply(y.counts, 1, sum)
    }
    
    y.prop = y.counts/apply(y.counts, 1, sum)  # compute proportions
    zdat = apply(y.prop, 2, glm.binomial.wrapper, g = g, W = W, weights = offset, 
        epsilon = 1e-06)  #estimate effect sizes and standard errors
    # zdat.ash = ash(zdat[3,],zdat[4,],df=2,method='fdr') #shrink the estimated
    # effects return(list(zdat=zdat,zdat.ash=zdat.ash))
    return(zdat)
}

glm.binomial.wrapper = function(y, g, W = NULL, ...) {
    if (is.null(W)) {
        y.glm = safe.quasibinomial.glm(y ~ g, ...)
    } else {
        y.glm = safe.quasibinomial.glm(y ~ g + W, ...)
    }
    return(c(get.coeff(y.glm), summary(y.glm)$dispersion))
}


# fill NAs with 0s (or other value)
fill.nas = function(x, t = 0) {
    x[is.na(x)] = t
    return(x)
}

# get estimates and standard errors from a glm object return NAs if not converged
get.coeff = function(x.glm) {
    c = as.vector(t(summary(x.glm)$coeff[, 1:2]))
    if (x.glm$conv) {
        return(c)
    } else {
        return(rep(NA, length(c)))
    }
}

# use glm to fit quasibinomial, but don't allow for underdispersion!
safe.quasibinomial.glm = function(formula, forcebin = FALSE, ...) {
    if (forcebin) {
        fm = glm(formula, family = binomial, ...)
    } else {
        fm = glm(formula, family = quasibinomial, ...)
        if (is.na(summary(fm)$dispersion) | summary(fm)$dispersion < 1) {
            fm = glm(formula, family = binomial, ...)
        }
    }
    return(fm)
}

# randomly subsample data for each gene gene: a vector of reads for one gene
# Nsamp: # of samples wanted
sampleingene = function(gene, Nsamp) {
    sample = sample(length(gene), Nsamp)
    return(c(gene[sample], sample))
}

# Randomly select samples counts: full count matrix Nsamp: # of samples wanted
# breaksample: flag, if select different samples for each gene
selectsample = function(counts, Nsamp, breaksample) {
    if (breaksample == FALSE) {
        subsample = sample(1:dim(counts)[2], Nsamp)
        counts = counts[, subsample]
        subsample = t(matrix(rep(subsample, dim(counts)[1]), ncol = dim(counts)[1]))
    } else {
        temp = t(apply(counts, 1, sampleingene, Nsamp = Nsamp))
        counts = temp[, 1:Nsamp]
        subsample = temp[, (Nsamp + 1):(2 * Nsamp)]
    }
    return(list(counts = counts, subsample = subsample))
}

# Use RUV to estimate confounding factor
RUV_factor = function(counts, args, null) {
    seq = newSeqExpressionSet(as.matrix(counts[as.logical(null), ]))
    if (sum(null) > 0) {
        controls = rownames(seq)
        differences = matrix(data = c(1:args$Nsamp, (args$Nsamp + 1):(2 * args$Nsamp)), 
            byrow = TRUE, nrow = 2)
        # seqRUV = RUVs(seq, controls, k=args$RUV.k, differences)
        seqRUV = RUVg(seq, controls, k = args$RUV.k)
        return(W = as.matrix(pData(seqRUV)))
    } else {
        return(W = NULL)
    }
}

# Use SVA to estimate confounding factor
SVA_factor = function(counts, condition, args, null) {
    mod1 = model.matrix(~condition)
    mod0 = cbind(mod1[, 1])
    
    if (args$nullpi > 0) {
        svseq = svaseq(counts, mod1, mod0, control = null)
    } else {
        svseq = svaseq(counts, mod1, mod0)
    }
    
    if (svseq$n.sv > 0) {
        return(W = svseq$sv)
    } else {
        return(W = NULL)
    }
}

# Get sebetahat from edgeR.glm (infer from betahat & pval)
edgeR_glmest = function(counts, condition, args) {
    design = model.matrix(~condition)
    y = DGEList(counts = counts + args$pseudocounts, group = condition)
    y = calcNormFactors(y)
    y = estimateGLMCommonDisp(y, design)
    y = estimateGLMTrendedDisp(y, design)
    y = estimateGLMTagwiseDisp(y, design)
    fit = glmFit(y, design)
    lrt = glmLRT(fit, coef = 2)
    betahat.edgeRglm = fit$coef[, 2]
    df.edgeRglm = length(condition) - 2
    tscore = qt(1 - lrt$table$PValue/2, df = df.edgeRglm)
    sebetahat.edgeRglm = abs(betahat.edgeRglm/tscore)
    
    return(list(betahat = betahat.edgeRglm, sebetahat = sebetahat.edgeRglm, df = df.edgeRglm))
}

# Get sebetahat from DESeq.glm (infer from betahat & pval)
DESeq_glmest = function(counts, condition, args) {
    cds = newCountDataSet(counts + args$pseudocounts, condition)
    cds = estimateSizeFactors(cds)
    cds = try(estimateDispersions(cds), silent = TRUE)
    if (class(cds) == "try-error") {
        betahat.DESeqglm = NA
        sebetahat.DESeqglm = NA
        df.DESeqglm = length(condition) - 2
    } else {
        fit1 = fitNbinomGLMs(cds, count ~ condition)
        fit0 = fitNbinomGLMs(cds, count ~ 1)
        betahat.DESeqglm = fit1[, 2]
        df.DESeqglm = length(condition) - 2
        tscore = qt(1 - nbinomGLMTest(fit1, fit0)/2, df = df.DESeqglm)
        sebetahat.DESeqglm = abs(betahat.DESeqglm/tscore)
    }
    return(list(betahat = betahat.DESeqglm, sebetahat = sebetahat.DESeqglm, df = df.DESeqglm))
}

# Get sebetahat from DESeq2 (infer from betahat & pval)
DESeq2_glmest = function(counts, condition, args) {
    cond = condition
    dds = DESeqDataSetFromMatrix(counts + args$pseudocounts, DataFrame(cond), ~cond)
    dds = estimateSizeFactors(dds)
    dds = estimateDispersions(dds)
    dds = nbinomWaldTest(dds)
    res = results(dds, cooksCutoff = FALSE)
    pvalue = res$pvalue
    betahat.DESeq2 = res$log2FoldChange
    df.DESeq2 = length(condition) - 2
    tscore = qt(1 - pvalue/2, df = df.DESeq2)
    sebetahat.DESeq2 = abs(betahat.DESeq2/tscore)
    
    return(list(betahat = betahat.DESeq2, sebetahat = sebetahat.DESeq2, df = df.DESeq2))
}


# Extract dataset for a specific tissue from the GTEx reads txt file Note: use
# GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct_new.txt, which removes
# the first 3 lines from
# GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.txt
readtissue = function(path, tissue) {
    tis = read.table(paste0(path, "/data/sample_tissue.txt"), header = TRUE)
    
    tissue.idx = grep(tissue, tis[, 2], fixed = TRUE)
    cols = rep("NULL", dim(tis)[1])
    cols[tissue.idx] = "numeric"
    
    data = read.table(paste0(path, "/gtex/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct_new.txt"), 
        colClasses = cols, header = TRUE)
    return(data)
} 


#' Generate "count" data that satisfies assumptions of SUCCOTASH.
#'
#' log_2(counts + 1) are normal with the factor-augmented regression model as the mean.
#'
#' @param n An integer. The number of rows
#' @param p An integer. The number of columns
#' @param k An integer. The number of confounders.
#' @param sig_daig A vector of positive numerics of length \code{p}. The variances of the columns.
#' @param sig_alpha A postive numeric. The standard devation of the confounder coefficients.
#' @param beta0 A positive numeric. The intercept for all regressions.
#' @param tau_seq A vector of positive numerics. The mixture standard deviations of beta
#' @param pi_seq A vector of positive numerics that sum to 1, pi_seqk[k]  is the
#'        probability of beta being drawn from a normal with variance tau_seq[k].
#' @param get_null A logical. Should beta be set to zero?
#' @return \code{counts}: An n by p matrix of counts
#'
#' \code{X}: An n by 2 design matrix.
#'
#' \code{beta}: A 2 by p matrix of coefficients.
#'
#' \code{Z}: An n by k matrix of hidden covariates.
#'
#' \code{alpha}: A k by p matrix of coefficients to the hidden confounders.
#'
#' \code{E}: An n by p matrix of error terms.
#'
#' \code{Y}: An n by p matrix of responses that are truly following the SUCCOTASH model.
theo_count_dat <- function(n, p, k, sig_diag, sig_alpha, beta0 = 10, get_null = TRUE, tau_seq = NULL, pi_seq = NULL)
{
    if(!get_null) {
        if(is.null(tau_seq) | is.null(pi_seq)) {
            warning("Need to specify tau_seq and pi_seq for nonnull data")
        }
        treat <- rep(1, length = n)
        treat[1:round(n/2)] <- 0
        X <- cbind(rep(1, n), treat)
        beta <- rbind(rep(beta0, p), draw_beta(pi_vals = pi_seq, tau_seq = tau_seq, p = p))
        Z <- matrix(rnorm(n * k), nrow = n, ncol = k)
        alpha <- matrix(rnorm(k * p, sd = sig_alpha), nrow = k, ncol = p)
        E <- matrix(rnorm(n * p), nrow = n, ncol = p) %*%  diag(sqrt(sig_diag))
        Y <- X %*% beta + Z %*% alpha + E
        counts <- round(2^Y - 1)
        counts[counts < 0] <- 0
    }
    else if(get_null) {
        treat <- rep(1, length = n)
        treat[1:round(n/2)] <- 0
        X <- cbind(rep(1, n), treat)
        beta <- rbind(rep(beta0, p), rep(0, length = p))
        Z <- matrix(rnorm(n * k), nrow = n, ncol = k)
        alpha <- matrix(rnorm(k * p, sd = sig_alpha), nrow = k, ncol = p)
        E <- matrix(rnorm(n * p), nrow = n, ncol = p) %*%  diag(sqrt(sig_diag))
        Y <- X %*% beta + Z %*% alpha + E
        counts <- round(2^Y - 1)
        counts[counts < 0] <- 0
    }
    
    return(list(counts = counts, X = X, beta = beta, Z = Z, alpha = alpha, E = E, Y = Y))
}


#' Generate "count" data that satisfies assumptions of SUCCOTASH.
#'
#' Same as theo_count_dat but now the log-counts-per-million are normal with the factor-augmented regression model as the mean.
#' Not the log_2(counts + 1) as before.
#'
#' @param n An integer. The number of rows
#' @param p An integer. The number of columns
#' @param k An integer. The number of confounders.
#' @param sig_daig A vector of positive numerics of length \code{p}. The variances of the columns.
#' @param sig_alpha A postive numeric. The standard devation of the confounder coefficients.
#' @param beta0 A positive numeric. The intercept for all regressions.
#' @param tau_seq A vector of positive numerics. The mixture standard deviations of beta
#' @param pi_seq A vector of positive numerics that sum to 1, pi_seqk[k]  is the
#'        probability of beta being drawn from a normal with variance tau_seq[k].
#' @param get_null A logical. Should beta be set to zero?
#' @param R_i a vector of length n. The total read count for individual i. Defaults to 10^7 for each individual.
#'
#' @return \code{counts}: An n by p matrix of counts
#'
#' \code{X}: An n by 2 design matrix.
#'
#' \code{beta}: A 2 by p matrix of coefficients.
#'
#' \code{Z}: An n by k matrix of hidden covariates.
#'
#' \code{alpha}: A k by p matrix of coefficients to the hidden confounders.
#'
#' \code{E}: An n by p matrix of error terms.
#'
#' \code{Y}: An n by p matrix of responses that are truly following the SUCCOTASH model.
theo_count_dat_cpm <- function(n, p, k, sig_diag, sig_alpha, beta0 = 10, get_null = TRUE, tau_seq = NULL, pi_seq = NULL, R_i = NULL)
{
    if(is.null(R_i)) {
        R_i <- 10^7
    }
    if(!get_null) {
        if(is.null(tau_seq) | is.null(pi_seq)) {
            warning("Need to specify tau_seq and pi_seq for nonnull data")
        }
        treat <- rep(1, length = n)
        treat[1:round(n/2)] <- 0
        X <- cbind(rep(1, n), treat)
        beta <- rbind(rep(beta0, p), draw_beta(pi_vals = pi_seq, tau_seq = tau_seq, p = p))
        Z <- matrix(rnorm(n * k), nrow = n, ncol = k)
        alpha <- matrix(rnorm(k * p, sd = sig_alpha), nrow = k, ncol = p)
        E <- matrix(rnorm(n * p), nrow = n, ncol = p) %*%  diag(sqrt(sig_diag))
        Y <- X %*% beta + Z %*% alpha + E
        counts <- 2^Y - 0.5
        counts <- round((R_i / rowSums(counts)) * counts)
        counts[counts < 0] <- 0
    }
    else if(get_null) {
        treat <- rep(1, length = n)
        treat[1:round(n/2)] <- 0
        X <- cbind(rep(1, n), treat)
        beta <- rbind(rep(beta0, p), rep(0, length = p))
        Z <- matrix(rnorm(n * k), nrow = n, ncol = k)
        alpha <- matrix(rnorm(k * p, sd = sig_alpha), nrow = k, ncol = p)
        E <- matrix(rnorm(n * p), nrow = n, ncol = p) %*%  diag(sqrt(sig_diag))
        Y <- X %*% beta + Z %*% alpha + E
        counts <- 2^Y - 0.5
        counts <- round((R_i / rowSums(counts)) * counts)
        counts[counts < 0] <- 0
    }
    
    return(list(counts = counts, X = X, beta = beta, Z = Z, alpha = alpha, E = E, Y = Y))
}



#args_val <- list()
#args_val$tissue <- "Lung"
#args_val$Nsamp <- 10
#args_val$path <- "../data/"
#args_val$Ngene <- 100

#args_val$pi_seq <- c(0.5, 0.3, 0.2)
#args_val$tau_seq <- c(0, 5, 10)
#args_val$get_null <- FALSE
#args_val$Nconfounders <- 5


#n <- dfargs$Nsamp
#p <- dfargs$Ngene
#k <- dfargs$Nconfounders
#sig_diag <- dfargs$sig_diag
#sig_alpha <- dfargs$sig_alpha
#beta0 <- dfargs$beta0
#get_null <- dfargs$get_null
#tau_seq <- dfargs$tau_seq
#pi_seq <- dfargs$pi_seq
