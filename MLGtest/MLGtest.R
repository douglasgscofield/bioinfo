# The MLGtest is a bootstrap test that examines whether the number of unique
# multilocus genotypes (MLGs) observed within a sample from a population is
# significantly different from the number expected based on observed allele
# frequencies.  The test compares the observed number of MLGs against a
# bootstrap distribution of number of MLGs for the given sample size generated
# using allele frequencies given for the observed population.  An arbitrary
# number of loci and alleles may be used for each population.
#
# This test is likely overconservative.

# allele.freqs is a data frame with four columns:
#
#   population: the population ID, identical for all entries for a population
#
#   locus:      locus ID, identical for all allele entries for a locus; may
#               be reused in different populations
#
#   allele:     allele ID, unique for each allele within a locus, may be
#               reused for different alleles
#
#   freq:       population frequency of the allele within the locus
#
#
# observed is a data frame with three columns:
#
#   population: the population ID, matching the population ID used in 
#               allele.freqs
#
#   n:          the number of individuals genotyped in the population
#
#   observed:   the number of unique multilocus genotypes observed in
#               the population


MLGtest <- function(allele.freqs=read.delim('allele-freqs.txt'), 
                    MLG.observed=read.delim('observed.txt'),
                    quantiles=c(0,0.01,0.05,0.1,0.5,0.9,0.95,0.99,1),
                    ploidy=2,
                    N.iter=1000,
                    check=TRUE,
                    include.dist=FALSE)
{
    if (check)
        allele.freqs <- MLG.check.allele.freqs(allele.freqs, drop.bad=TRUE, verbose=TRUE)
    cat("MLGtest:",nrow(MLG.observed),"populations,",N.iter,"Monte Carlo iterations\n")
    ans <- list()
    ans.df <- data.frame()
    for (pop in unique(allele.freqs$population)) {
        cat("Population",pop,"...\n")
        ans[[as.character(pop)]] <- list()
        pfreqs <- subset(allele.freqs, population == pop)
        pobs <- unlist(subset(MLG.observed, population == pop))
        dist <- pop.dist.genotypes(pfreqs, N.iter=N.iter, n=pobs['n'], ploidy=ploidy)
        ans.df <- rbind(ans.df, data.frame(pop=pop, n=n, ploidy=ploidy, N.iter=N.iter,
                                           MLG.mean=mean(dist), MLG.var=var(dist),
                                           MLG.observed=pobs['observed'],
                                           P.observed=sum(pobs['observed'] >= dist) / N.iter))
        ans[[pop]][['population']] <- pop
        ans[[pop]][['ploidy']] <- ploidy
        ans[[pop]][['n']] <- pobs['n']
        ans[[pop]][['N.iter']] <- N.iter
        ans[[pop]][['MLG.sim.mean']] <- mean(dist)
        ans[[pop]][['MLG.sim.var']] <- var(dist)
        qqq <- quantile(dist, quantiles)
        names(qqq) <- as.character(quantiles)
        ans[[pop]][['MLG.sim.quantiles']] <- qqq
        ans[[pop]][['MLG.observed']] <- pobs['observed']
        ans[[pop]][['P.observed']] <- (sum(pobs['observed'] >= dist) / N.iter)
        if (include.dist)
            ans[[pop]][['dist']] <- dist
    }
    for (pop in names(ans)) {
        cat(paste('=== Population', pop, '===\n'))
        for (v in c('ploidy','n','MLG.observed','N.iter','MLG.sim.mean','MLG.sim.var','P.observed'))
            cat(paste(v, '=', ans[[pop]][[v]], '\n'))
        qqq <- ans[[pop]][['MLG.sim.quantiles']]
        cat('MLG.sim.quantiles :\n')
        for (i in 1:length(qqq))
            cat(paste('   ', names(qqq)[i],'=', qqq[i], '\n'))
        cat('\n')
    }
    invisible(list(ans=ans, ans.df=ans.df))
}

# expect a single population

pop.dist.genotypes <- function(pdat, N.iter=1, n, ploidy=2)
{
    stopifnot(length(unique(pdat$population)) == 1) 
    stopifnot(! missing(n))
    summlist <- pop.genotype.summary(pdat)
    #if ("warn.conflicts" %in% names(as.list(args(attach))))
    #    attach(summlist, warn.conflicts=FALSE)
    #else attach(summlist)
    ans <- rep(0, N.iter)
    gt.allele.index <- matrix(rep(0, ploidy*summlist$nloci), ploidy, summlist$nloci)
    genotypes <- rep("", n)
    gt <- rep("", ploidy)

    n.uniform <- N.iter * n * ploidy * summlist$nloci
    uniform.deviate <- runif(n.uniform)

    for (trial in 1:N.iter) {
        u.index.2 <- ((trial - 1) * (n * ploidy * summlist$nloci))
        for (i in 1:n) {
            # step through loci, choosing a genotype
            u.index.1 <- u.index.2 + ((i - 1) * (ploidy * summlist$nloci))
            for (hap in 1:ploidy) {
                u.index.0 <- u.index.1 + ((hap - 1) * summlist$nloci)
                for (loc in 1:summlist$nloci) {
                    u.index <- u.index.0 + loc
                    u <- uniform.deviate[u.index]
                    #uniform.fresh[u.index] <- uniform.fresh[u.index] - 1
                    for (al in 1:summlist$nalleles[loc]) {
                        if (u <= summlist$cumfreq[[loc]][al]) {
                            gt.allele.index[hap,loc] <- al
                            break
                        }
                    }
                }
            }
            # we don't care about the order, so sort 'em
            gt.allele.index <- apply(gt.allele.index, 2, sort)
            for (hap in 1:ploidy) {
                gt[hap] <- paste(collapse="", gt.allele.index[hap,])
            }
            genotypes[i] <- paste(collapse="_", gt)
        }
        ans[trial] <- length(unique(genotypes))
    }
    ans
} 
pop.genotype.summary <- function(pdat)
{
    if (length(unique(pdat$population)) > 1) 
        stop("pop.genotype.summary needs just one population")
    loci <- unique(as.character(pdat$locus))
    nloci <- length(loci)
    alleles <- list()
    nalleles <- rep(0,nloci)
    cumfreq <- list()
    for (lo in 1:nloci) {
        spd <- subset(pdat, locus == lo)
        alleles[[lo]] <- as.character(spd$allele)
        if (any(duplicated(alleles[[lo]])))
            stop("nonunique alleles found")
        nalleles[lo] <- length(alleles[[lo]])
        cumfreq[[lo]] <- cumsum(spd$freq)
    }
    list(nloci=nloci, loci=loci, nalleles=nalleles, alleles=alleles, cumfreq=cumfreq)
}

MLG.check.allele.freqs <- function(dat, drop.bad = TRUE, verbose=TRUE) 
{
    newdat <- check.locus.sums(dat, verbose=verbose)
    if (drop.bad) {
        if (nrow(newdat) < nrow(dat) && verbose) 
            cat(paste("there were locus frequences that didn't sum to 1\n"))
        return(newdat)
    }
    else return(dat)
}

check.locus.sums <- function(dat, epsilon=1e-5, verbose=TRUE)
{
    newdat <- dat
    for (pop in unique(dat$population)) {
        sd <- subset(dat, population == pop)
        for (loc in unique(sd$locus)) {
            lsd <- subset(sd, locus == loc)
            sumlsd <- sum(lsd$freq)
            if (abs(1 - sumlsd) > epsilon) {
                if (verbose) {
                    cat(paste(lsd$freq, "\n"))
                    cat(paste("1 - sumlsd =", 1 - sumlsd, "\n"))
                    cat(paste("pop", pop, "locus", loc, "sumlsd =", sumlsd, "\n"))
                }
                # drop the locus
                newdat <- subset(newdat, population != pop | locus != loc)
            }
        }
    }
    return(newdat)
}

