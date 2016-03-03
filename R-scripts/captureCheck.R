# Compare chloroplast SNP frequencies between pools from different plant
# species and check for the possibility of chloroplast capture by one or
# more individuals.  The method is based on frequencies from pooled data so
# probably of limited usefulness for individuals.

# Copyright (c) 2013 Douglas G. Scofield, douglas.scofield@ebc.uu.se
# Uppsala University (Umea University at the time of authorship)

captureCheck = function(snps,     # data.frame of SNPs for a group of species
                         T = 0.95, # threshold frequency for fixed SNP
                         Tp = 0.1, # threshold low frequency for contaminant SNP
                         alt.spp = "T", # spp to use alternate fixation threshold
                         alt.T = 0.90, # alternate fixation threshold
                         output = c("report", "table"))
{
    output = match.arg(output)
    on.exit(options(stringsAsFactors=FALSE))
    tab = data.frame(Ref=I(""), 
                     Contam=I(""),
                     OtherSpp=I(""),
                     fixed.C=0,
                     fixed.C.in.R=0,
                     fixed.C.in.R.lowfreq=0,
                     fixed.C.in.R.medianfreq=0,
                     fixed.CO=0,
                     fixed.CO.in.R=0,
                     fixed.CO.in.R.lowfreq=0,
                     fixed.CO.in.R.medianfreq=0)[0, ]
    snps$fixed = snps$SNPAllele != snps$ReferenceAllele & snps$Proportion >= T
    if (alt.spp != "") {
        # An alternate spp was specified with its own fixation threshhold
        # Important to consider because contamination screws up "fixation"
        alt.fixed = snps$species %in% alt.spp & snps$SNPAllele != snps$ReferenceAllele & snps$Proportion >= alt.T
        snps$fixed = snps$fixed | alt.fixed
    }
    spp = sort(unique(snps$species))
    locs.fixed = list()
    nRef = function(x) paste(sep=".", "Ref", x)
    nContamFixed = function(x) paste(sep=".", "Contam", x, "fixed")
    nContam = function(x) paste(sep=".", "Contam", x)
    nlowfreq = function(x) paste(sep=".", x, "lowfreq")
    for (R in spp) {  # is R contaminated?
        cat("Checking for contamination in", R)
        # find set of SNPs in R
        snps.R.table = subset(snps, species == R & Proportion > 0)
        rownames(snps.R.table) = snps.R.table$Position
        snps.R.fixed = subset(snps.R.table, fixed)$Position
        snps.R.lowfreq = subset(snps.R.table, Proportion <= Tp)$Position
        snps.R = snps.R.table$Position
        cat("  (", length(snps.R), "SNPs,", 
            length(snps.R.fixed), "fixed >=", if (R %in% alt.spp) alt.T else T,
            length(snps.R.lowfreq), "low <=", Tp, 
            ")\n")
        locs.fixed[[nRef(R)]] = list()
        for (C in setdiff(spp, R)) {  # here, C is the potential contaminant
            cat("-- Is", C, "a contaminant?\n")
            O = setdiff(spp, c(R, C))  # species other than C and R
            # produce set of fixed SNPs
            fixed = with(subset(snps, fixed), table(species, Position))
            # find set of SNPs that are fixed in C that are not fixed in R or O
            fixed.C = fixed[C, ] & !apply(fixed[c(R, O), ], 2, any)
            fixed.C = as.numeric(colnames(fixed[, fixed.C]))
            locs.fixed[[nRef(R)]][[nContamFixed(C)]] = fixed.C
            cat("----", length(fixed.C), "SNPs fixed only in", C, "\n")
            fixed.C.in.R = fixed.C %in% snps.R
            locs.fixed[[nRef(R)]][[nContam(C)]] = fixed.C[fixed.C.in.R]
            cat("------", sum(fixed.C.in.R), "of these (",
                round(sum(fixed.C.in.R)*100/length(fixed.C), 2), 
                "%) also found in", R, "\n")
            fixed.C.in.R.lowfreq = fixed.C %in% snps.R.lowfreq
            locs.fixed[[nRef(R)]][[nlowfreq(nContam(C))]] = fixed.C[fixed.C.in.R.lowfreq]
            cat("------", sum(fixed.C.in.R.lowfreq), "of these (",
                round(sum(fixed.C.in.R.lowfreq)*100/length(fixed.C), 2), 
                "%) are of low frequency ( <= ", Tp, ") in", R, "\n")
            median.freq.C = if (any(fixed.C.in.R))
                median(snps.R.table[as.character(fixed.C[fixed.C.in.R]), ]$Proportion)
            else
                0

            # find set of SNPs that are fixed in C and O but not R
            CO = paste(collapse=",", c(C, O))
            fixed.CO = apply(fixed[c(C, O), ], 2, all) & !fixed[R, ]
            fixed.CO = as.numeric(colnames(fixed[, fixed.CO]))
            locs.fixed[[nRef(R)]][[nContamFixed(CO)]] = fixed.CO
            cat("----", length(fixed.CO), "SNPs fixed in", CO, "and not", R, "\n")
            fixed.CO.in.R = fixed.CO %in% snps.R
            locs.fixed[[nRef(R)]][[nContam(CO)]] = fixed.CO[fixed.CO.in.R]
            cat("------", sum(fixed.CO.in.R), "of these (",
                round(sum(fixed.CO.in.R)*100/length(fixed.CO), 2), 
                "%) also found in", R, "\n")
            fixed.CO.in.R.lowfreq = fixed.CO %in% snps.R.lowfreq
            locs.fixed[[nRef(R)]][[nlowfreq(nContam(CO))]] = fixed.CO[fixed.CO.in.R.lowfreq]
            cat("------", sum(fixed.CO.in.R.lowfreq), "of these (",
                round(sum(fixed.CO.in.R.lowfreq)*100/length(fixed.CO), 2), 
                "%) are of low frequency ( <=", Tp, ") in", R, "\n")
            median.freq.CO = if (any(fixed.CO.in.R))
                median(snps.R.table[as.character(fixed.CO[fixed.CO.in.R]), ]$Proportion)
            else
                0

            l = list(Ref=I(R),
                     Contam=I(C),
                     OtherSpp=I(paste(collapse=",", O)),
                     fixed.C=length(fixed.C),
                     fixed.C.in.R=sum(fixed.C.in.R),
                     fixed.C.in.R.lowfreq=sum(fixed.C.in.R.lowfreq),
                     fixed.C.in.R.medianfreq=median.freq.C,
                     fixed.CO=length(fixed.CO),
                     fixed.CO.in.R=sum(fixed.CO.in.R),
                     fixed.CO.in.R.lowfreq=sum(fixed.CO.in.R.lowfreq),
                     fixed.CO.in.R.medianfreq=median.freq.CO)
            tab = rbind(tab, l)
        }
    }
    # add locations of fixed Contam and Contam+OtherSpp SNPs
    attr(tab, "locs.fixed") = locs.fixed
    ####
    tab
}

