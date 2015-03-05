###########################################################################
#
# copyright (c) 2011 Douglas G. Scofield, douglasgscofield@gmail.com
#
# Use as you see fit, just include the above attribution.  No warranty is
# implied or assumed with this code.


###########################################################################
# convertHapSNPsToArlequin
#
# Convert a file containing population-specific *haplotype SNPs* in a simple 
# format to Arlequin project file format.  Input is a single file with or
# without a header line (depending upon option has.header) containing columns
# in the following order: 
#
#    samplename popname loc1 loc2 loc3 ... locN
#
# Genotypes at each locus should be SNPs specified using DNA codes understood
# by Arlequin.  The code for missing data is supplied via the option
# missing.data.  Unique haplotypes are determined by concatenating genotypes
# together with no space between loci; Arlequin currently (version 3.5)
# considers each nucleotide to be a separate locus, which is appropriate for
# SNPs.  Haplotypes are named "hap<N>" where N is incremented as unique
# haplotypes are encountered.
#
# Output conforms to Arlequin's requirements for a project file; see Arlequin's
# documentation for more information.  Output includes [Profile] options that
# are specific to DNA haplotypes with concatenated genotypes:  
#
#   DataType=DNA
#   GenotypicData=0
#   LocusSeparator=NONE
#
# Individual samples are grouped within populations (which Arlequin calls
# samples).
#
# Note that this tool only formats genotypic data and can provide information
# on very limited types of population structure; Arlequin's project file format
# can contain much more information and a wider variety of data input formats
# (e.g., distance matrices) than the genotypic format produced here.
#
# The website for the current version of Arlequin (3.5) may be found at
#
#   http://popgen.unibe.ch/software/arlequin35/
#
#
# Arguments:
#
# infile            Name of input file containing sample records, with each
#                   line containing *in order* tab-separated columns for the 
#                   sample name, population for the sample, and genotypes at 
#                   any number of haploid loci.  If the first line of the file 
#                   has column headings, then use has.header=TRUE.
#                   (default "example_data.txt")
#
# outfile           Name of output file as an Arlequin project file.  If the
#                   outfile name is empty "", then output is produced within
#                   the R console and no file is created.  (default "<infile>.arp").
#
# has.header        If TRUE, the first line of infile is column headings.
#                   (default FALSE)
#
# missing.data      The code used within the locus columns to indicate missing
#                   data.  This does not result in any special processing of
#                   the locus values, rather it is directly output into the
#                   the Arlequin project file in the [Profile] option 
#                   MissingData.  (default "?")
#
# title             Title within the Arlequin project file, output in the 
#                   [Profile] option Title.  (default "<infile> haplotypes")
#
# sample.method     If "by.haplotype" (the default), then unique sample
#                   haplotypes within each population are grouped and given a
#                   number of occurrences, e.g., if haplotype "ACGT" is given 
#                   the name "hap1" and appears twice within population "pop1", 
#                   then the output for the population is
#
#                     SampleName="pop1"
#                     SampleSize=2
#                     SampleData={
#                       hap1	2
#                     }
#
#                   If "by.individual", then each individual within each
#                   population is output, using the sample identifier:
#
#                     SampleName="pop1"
#                     SampleSize=2
#                     SampleData={
#                       samp01	1
#                       samp02	1
#                     }
#
#                   In general you may want to stick with "by.haplotype".  If
#                   you choose "by.individual", then your sample names must not
#                   have whitespace within them, and there are additional
#                   considerations required for your Arlequin analyses to
#                   recognize different individuals as having haplotypes that
#                   are effectively the same  See the Arlequin documentation
#                   for more information.
#
# haplotype.method  If "section" (the default), then haplotype genotypes are 
#                   output in the [[HaplotypeDefinition]] section of the 
#                   Arlequin project file:
#
#                     [[HaplotypeDefinition]]
#                       HaplListName="<infile> haplotype list"
#                       HaplList={
#                         hap1	ACGT
#                       }
#
#                   If "file", then a separate haplotype file is created
#                   named "<infile>.hap", containing one line with two
#                   tab-separated columns for each haplotype and no column 
#                   headings:
#
#                         hap1	ACGT
#                         ....
#
#                   and this file is referred to within the 
#                   [[HaplotypeDefinition]] section of the project file:
#
#                     [[HaplotypeDefinition]]
#                       HaplListName="<infile> haplotype list"
#                       HaplList = EXTERN "<infile>.hap"
#
#                   This file must be available to Arlequin when it reads
#                   the project file, in the same folder as the project file.
#
#                   If "inline", then haplotypes are output with each sample
#                   line, and there is no separate haplotype section or file:
#
#                     SampleData={
#                       hap1	2  ACGT
#                     }
#
# structure.method  If "none" (the default), then no [[Structure]] section
#                   is output in the Arlequin project file.
#
#                   If "by.population", then a [[Structure]] section is
#                   produced with each population forming its own group:
#
#                     [[Structure]]
#                       StructureName="<infile> structure = by.population"
#                       NbGroups=2
#                       Group= {
#                         "pop1"
#                       }
#                       Group= {
#                         "pop2"
#                       }
#
#                   If "all", then a [[Structure]] section is produced 
#                   with all populations lumped into a single group:
#
#                     [[Structure]]
#                       StructureName="<infile> structure = all"
#                       NbGroups=1
#                       Group= {
#                         "pop1"
#                         "pop2"
#                       }
#


convertHapSNPsToArlequin = function(
                infile="example_data.txt", 
                outfile=paste(sep="", infile, ".arp"),
                has.header=FALSE,
                missing.data="?",
                title=paste(infile, "haplotypes"),
                sample.method=c("by.haplotype", "by.individual"),
                haplotype.method=c("section", "file", "inline"),
                structure.method=c("none", "by.population", "all")
                )
{
  .version = "convertHapSNPsToArlequin v. 0.1"

  on.exit(suppressWarnings(sink()))  # try to clean up nicely if we terminate early
  locus.separator = ""  # must be empty string, change output of LocusSeparator= if not

  .qu = function(x)
  { # place double quotes around a string
    paste(sep="", collapse="", "\"", x, "\"")
  }
  .arlcat = function(indent, ...) 
  { # print at a given indent level, with no separator
    indent.string = "  "
    while (indent > 0) {
      cat(sep="", indent.string)
      indent = indent - 1
    }
    cat(sep="", ...)
  }

  sample.method = match.arg(sample.method)
  haplotype.method = match.arg(haplotype.method)
  structure.method = match.arg(structure.method)

  if (class(infile) == "character")
    dat = read.delim(infile, header=has.header, sep="\t", as.is=TRUE)
  else 
    stop("don't recognize filename")

  names(dat)[1:2] = c("sample", "pop")

  # form haplotype of each, tag with haplotype name
	pops = sort(unique(dat$pop))
  Ncol = ncol(dat)
	dat$hap = apply(dat[,3:Ncol], 1, function(x) paste(collapse=locus.separator, x))
	uhaps = unique(dat$hap)
	names(uhaps) = paste(sep="", "hap", 1:length(uhaps))
	dat$hapnames = sapply(dat$hap, function(x) paste(sep="", "hap", which(x == uhaps)))
  NbSamples = length(pops)

	haps.df = data.frame(name=I(names(uhaps)), haplotype=I(uhaps))
  rownames(haps.df) = haps.df$name

	haplotype.file = paste(sep="", infile, ".hap")
  if (haplotype.method == "file") {
    write.table(file=haplotype.file, quote = FALSE, haps.df, col.names=FALSE, 
                row.names=FALSE, sep="\t");
  }

  if (outfile != "") # redirect output to outfile
    sink(file=outfile)

  .arlcat(0, "# Produced using ", .version, "\n")
  .arlcat(0, "[Profile]\n")
  .arlcat(2, "Title=", .qu(title), "\n")
  .arlcat(2, "NbSamples=", length(pops), "\n")
  .arlcat(2, "DataType=DNA\n")
  .arlcat(2, "GenotypicData=0\n")
  .arlcat(2, "LocusSeparator=NONE\n")
  .arlcat(2, "MissingData='", missing.data, "'\n")
  # .arlcat(2, "CompDistMatrix=1\n")  # double check

  .arlcat(0, "[Data]\n")

  if (haplotype.method %in% c("file", "section")) {
    haplotype.list.name = paste(infile, "haplotype list")
    .arlcat(1, "[[HaplotypeDefinition]]\n")
    .arlcat(2, "HaplListName=", .qu(haplotype.list.name), "\n")
    if (haplotype.method == "file") {
      .arlcat(2, "HaplList = EXTERN ", .qu(haplotype.file), "\n")
    } else if (haplotype.method == "section") {
      .arlcat(2, "HaplList={\n")
      for (r in 1:nrow(haps.df)) {
        .arlcat(3, haps.df[r, "name"], "\t", haps.df[r, "haplotype"], "\n")
      }
      .arlcat(2, "}\n")
    }
  }

  .arlcat(1, "[[Samples]]\n")

	for (p in pops) {

		.arlcat(2, "SampleName=",.qu(p),"\n")
		sdat = subset(dat, pop == p)
		.arlcat(2, "SampleSize=", nrow(sdat), "\n")
		.arlcat(2, "SampleData={\n")

    if (sample.method == "by.haplotype") {

      haps = table(sdat$hapnames)
      for (nm in names(haps)) {
        if (haplotype.method == "inline") {
          .arlcat(3, nm, "\t", haps[nm], "\t", haps.df[nm, "haplotype"], "\n")
        } else {
          .arlcat(3, nm, "\t", haps[nm], "\n")
        }
      }
    } else if (sample.method == "by.individual") {
      for (r in 1:nrow(sdat)) {
        if (haplotype.method == "inline") {
          .arlcat(3, sdat[r, "sample"], "\t", 1, "\t", sdat[r, "hap"], "\n")
        } else {
          .arlcat(3, sdat[r, "sample"], "\t", 1, "\n")
        }
      }
    }
		.arlcat(2, "}\n")
  }

  if (structure.method != "none") {
    .arlcat(1, "[[Structure]]\n")
    structure.name = paste(infile, "structure =", structure.method)
    .arlcat(2, "StructureName=", .qu(structure.name), "\n")
    if (structure.method == "by.population") {
      NbGroups = length(pops)
      .arlcat(2, "NbGroups=", NbGroups, "\n")
      for (p in pops) {
        .arlcat(2, "Group= {\n")
        .arlcat(3, .qu(p), "\n")
        .arlcat(2, "}\n")
      }
    } else if (structure.method == "all") {
      NbGroups = 1
      .arlcat(2, "NbGroups=", NbGroups, "\n")
      .arlcat(2, "Group= {\n")
      for (p in pops) {
        .arlcat(3, .qu(p), "\n")
      }
      .arlcat(2, "}\n")
    }
	}

  if (outfile != "") # stop redirecting output
    sink()

}


