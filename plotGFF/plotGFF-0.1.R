# (c) Douglas G. Scofield, 2011, McGill University
# douglasgscofield@gmail.com
# Use/modify/redistribute freely, keeping the above attribution
#
# plotGFF: graphic plot of GFF annotation
#
# An extended version built upon the already very useful plotRanges() available
# at
# http://www.bioconductor.org/help/course-materials/2009/SeattleNov09/IRanges/IRangesOverview.R
#
# Uses the following BioConductor packages.  Install BioConductor packages with
# biocLite("") after running source("http://bioconductor.org/biocLite.R")

library(rtracklayer)  # getting this with biocLite installs the other two...
library(IRanges)
library(Biostrings)

# x is an object containing the contents of a GFF file, created by import.gff()
# from the BioConductor package rtracklayer.

# xlim gives a range of coordinates from the GFF to which the plot is (more or
# less) restricted.

# sites is a list of sites (indicated by vertical lines) to add in tracks below
# the GFF-based annotation.  I use it to add SNP annotation (including optional
# SNP frequency) to a plot.  It can have one of a few formats.
#
#    sites=list(SNPs=c(5,10,15)) 
#       plot vertical lines at each coordinate, labeling the plot area "SNPs"
#
#    sites=list(SNPs=list(x=c(5,10,15), y=c(0.1,0.4,0.6)))
#       plot vertical lines at coordinates in x, with heights scaled 
#       following [0,1]-range values given in y.
#
#    sites=list(add.rect=TRUE, ...)
#       with either of the above, plot a rectangle around the set of vertical 
#       lines
#
# There can be multiple such lists wtihin the sites argument, resulting in
# multiple tracks of sites being plotted below the GFF.

# col.method indicates how the boxes should be colored: by type as indicated in
# the GFF, by strand orientation, or with a single color ("default")

# stack.method indicates how the boxes should be grouped horizontally and
# stacked vertically: by type from GFF annotation or by disjointBins analysis
# with Bioconductor.

# col is box color, recycled if necessary.  A titled key will be plotted if
# col.method is "type" or "strand".

# drop.source indicates to drop the 'source' type from the plotted annotation


plotGFF <-
function(x, 
         sites = NULL, # list(name1 = x1, name2 = x2,  ...)
         xlim = x, 
         main,
         col.method = c("type", "strand", "default"), 
         stack.method = c("type", "disjointBins"),
         col = c("black", "gray"), 
         sep = 0.2, 
         drop.source = TRUE,
         ...) 
{
  if (missing(main))
    main <- deparse(substitute(x))
  stack.method <- match.arg(stack.method)
  height <- 0.5
  sites.height <- 1.2
  site.height <- 0.9
  if (is(xlim, "RangedData")) {
    xlim <- c(min(start(xlim)), max(end(xlim)))
  }
  if (drop.source)  # drop source from types
    x <- subset(x, type != "source")
  col.method <- match.arg(col.method)
  typ <- unique(x$type)
  if (col.method == "default") {
    col <- col[1]
  } else if (col.method != "default" && ! all(col == c("black", "gray"))) {
    col <- col
  } else if (col.method == "strand") {
    col <- sapply(strand(x), 
                  function(.x) if (.x == "+") col[1] else 
                                 if (.x == "-") col[2] else "red")
  } else if (col.method == "type") {
    typecol <- rainbow(length(typ))
    names(typecol) <- typ
    col <- sapply(x$type, function(.x) typecol[as.character(.x)])
  }
  if (stack.method == "disjointBins") {
    bins <- disjointBins(IRanges(start(x), end(x) + 1))
  } else if (stack.method == "type") {
    bins <- match(x$type, typ)
  } else {
    stop("stack.method ", stack.method, " not implemented")
  }
  par(mar = c(4, 3, 1, 1), xpd = NA)
  plot.new()
  plot.window(xlim, c(0, max(bins + ifelse(is.null(sites), 0, 1))*(height + sep)))
  ybottom <- bins * (sep + height) - height
  if (! is.null(sites))
    ybottom <- ybottom + 1
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, 
       border = col, ...)
  if (! is.null(sites)) {
    # process options
    add.rect <- if (! is.null(sites$add.rect)) {
      t1 <- sites$add.rect
      sites$add.rect <- NULL
      t1
    } else FALSE
    site.y.range <- c(0, sites.height - sep)
    n.sites <- length(names(sites))
    st.height <- each.st.height <- (diff(site.y.range) / n.sites) * site.height
    st.y <- rev(seq(site.y.range[1], site.y.range[2] - each.st.height, length.out = n.sites))
    for (i in seq(along=names(sites))) {
      if (add.rect) rect(xlim[1], st.y[i], xlim[2], st.y[i] + each.st.height, 
                         lwd = 0.3, border = "blue")
      site.x <- if (is.list(sites[[i]])) sites[[i]]$x else sites[[i]]
      st.height <- if (is.list(sites[[i]]) && ! is.null(sites[[i]]$y))
        each.st.height * sites[[i]]$y
      else 
        each.st.height
      cat("site", names(sites)[i], "y0 = ", st.y[i], 
          " each.st.height =", each.st.height, 
          "max(st.height) =", max(st.height), "\n")
      segments(x0 = site.x, y0 = st.y[i], y1 = st.y[i] + st.height)
    }
    text(0, st.y + (each.st.height / 2), labels = names(sites), pos = 2)
  }
  if (col.method == "strand") {
    legend("topleft", title="Strand", pch=15, pt.cex=1.5, col=c(col,"red"), 
           legend=c("+","-","other"), bty="n")
  } else if (col.method == "type") {
    legend("topleft", title="Type", pch=c(rep(15, length(typ))), pt.cex=1.5, 
           col=typecol, legend=typ, bty="n")
  }
  title(main)
  axis(1)
}

