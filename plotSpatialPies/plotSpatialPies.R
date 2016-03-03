###########################################################################
#
# copyright (c) 2010-2011 Douglas G. Scofield, douglasgscofield@gmail.com
#
# Use as you see fit, include the above attribution.  No warranty is implied or
# assumed with this code.


###########################################################################
# plotSpatialPies
#
# Plot pie charts representing fractions of a code at particular locations.  Input is 
# one or two files containing named columns, and the plot can be controlled with a 
# number of options.  This was originally developed to plot proportions of particular
# genotypes at population locations, but is can now plot more general proportions of
# each 'code' at each site.
#
# Requires the library "plotrix" to plot floating pies.  At the R prompt, run
# install.packages("plotrix").
#
# As a quick intro, try plotSpatialPies("example_data.txt"), with the file containing
# lines representing samples at a collection of sites, and at least has columns named 
# (with any capitalization) 'site', 'code', 'lat', and 'long'.
#
# Arguments:
#
# infile            Name of file containing sample-by-sample records, with each line
#                   containing at least the columns 'site', 'code', 'lat' and 'long',
#                   with any capitalization.  For each unique site, the proportions of 
#                   the pie are determined from the proportions of each code 
#                   represented at the site.  The location of the pie for each site 
#                   is calculated from the mean lat and long among all samples for the 
#                   site.  'Save as' the file from an Excel sheet with appropriate 
#                   column names as "Tab-delimited text (.txt)".
#
# colors            Name of file containing at least the columns 'code' and 'color',
#                   with any capitalization, specifying the plotting color to be used 
#                   for each code in the infile.  'Save as' the file from an Excel 
#                   sheet with appropriate column names as "Tab-delimited text (.txt)".
#
# use.size          TRUE/FALSE, to scale the size of the pies according to sample size 
#                   at each site.  FALSE sets all pies to the same size. (default TRUE)
#
# use.latlong       TRUE/FALSE, to plot pies at calculated lat/long values.  If FALSE,
#                   an array of pies is plotted.  (default TRUE)
#
# reverse.long      TRUE/FALSE, to reverse longitude values from those given, e.g., 
#                   for Western Hemisphere. (default TRUE)
#
# add.legend        If not "", location of legend of code colors.  To click on the plot
#                   to set the location of the legend, specify the option as 
#                   add.legend=locator().  (default "topright")
#
# fill.method       If "table", use the colors file; if "rainbow", draw from a rainbow 
#                   palette.  If there is an error opening the colors file, then 
#                   "rainbow" is used.  (default "table")
#
# ...               Additional arguments passed to pie() or floating.pie()
#
#
# Change log:
#
# v. 0.1: initial release
# v. 0.2: more versatile avoidance of legend plot (add.legend=NA or =NULL)


plotSpatialPies <- function(infile="example_data.txt", 
                         colors="example_colors.txt",
                         use.size=TRUE,
                         use.latlong=TRUE,
                         reverse.long=TRUE,
                         add.legend="topright",
                         fill.method=c("table", "rainbow"),
                         ...)
{

  .version = "plotSpatialPies v. 0.2"

  if (! require(plotrix))
    stop("library plotrix not available, please run 'install.packages(\"plotrix\"'")

  method <- "pie"

  fill.method <- match.arg(fill.method)
  if (fill.method == "table") {
    code.colors = try(read.delim(colors, as.is=TRUE))
    if (class(code.colors) == "try-error") {
      cat(sep="", "*** Problem opening color file \"", colors, "\", switching to fill.method=\"rainbow\"\n")
      fill.method = "rainbow"
    }
  }

  dat = read.delim(infile, sep="\t", header=TRUE)
  names(dat) = tolower(names(dat))
  if (! all(c("site", "lat", "long", "code") %in% names(dat)))
    stop("infile must contain columns (ignoring capitalization) named 'site', 'lat', 'long', and 'code'\n")

  tab = with(dat, table(site, code))

  tab <- t(tab)

  n.samples <- apply(tab, 2, sum)

  fill.colors <- rep("", nrow(tab))
  names(fill.colors) <- rownames(tab)

  if (fill.method == "table") {

    names(code.colors) = tolower(names(code.colors))

    if (! all(c("code", "color") %in% names(code.colors)))
      stop("colors file must contain columns (ignoring capitalization) named 'code' and 'color'\n")

    codes.missing <- which(! code.colors$code %in% rownames(tab))
    if (length(codes.missing) > 0) 
      cat("*** Info: more colors provided than codes in the data:", code.colors$code[codes.missing])
    colors.missing <- which(! rownames(tab) %in% code.colors$code)
    if (length(colors.missing) > 0) {
      stop("some codes not given colors: ", rownames(tab)[colors.missing])
    }

    fill.colors[code.colors$code] <- code.colors$color

  } else if (fill.method == "rainbow") {

    fill.colors[] <- rainbow(length(fill.colors))

  } else {

    stop("non-table and rainbow fill methods not implemented yet")

  }

  # create site information

  t.site = with(dat, unlist(names(split(long, site))))
  t.long = with(dat, unlist(lapply(split(long, site), function(x) mean(x, na.rm=TRUE))))
  t.lat = with(dat, unlist(lapply(split(lat, site), function(x) mean(x, na.rm=TRUE))))
  sites = data.frame(site=t.site, long=if (reverse.long) -t.long else t.long, lat=t.lat)


  if (method == "bar") {

    stop("bar method not implemented yet")

  } else if (method == "pie") {

    n.pie <- ncol(tab)

    if (! use.latlong) {
      n.pie.col <- min(10, ceiling(n.pie / 2))
      n.pie.row <- ceiling(n.pie/n.pie.col)
      par(mfrow=c(n.pie.row, n.pie.col), mar=c(1,1,1,1), bty="n", xpd=NA)
    } else {
      par(xpd=NA)
      plot(sites$long, sites$lat, pch='.', col="white", bty="n", xlab="Longitude", ylab="Latitude")
    }
    for (p in 1:n.pie) {
      xvals <- tab[,p]
      xvals <- xvals[xvals > 0]
      colvals <- fill.colors[names(xvals)]
      N <- sum(xvals)
      # scale pie radius according to plot range, if necessary; this will need fiddling
      #max.radius = 0.95
      max.radius = min(c(0.95, diff(range(sites$lat)/10), diff(range(sites$long)/10)))
      radius <- if (use.size) {
        # 0.25 + (0.7 * N / max(n.samples))
        (max.radius*0.25) + ((max.radius*0.75) * sqrt(N) / sqrt(max(n.samples)))
        # 0.95 * sqrt(N) / sqrt(max(n.samples))
      } else {
        max.radius
      }

      fill.args <- list(col=colvals, border=NA)
      site <- colnames(tab)[p]
      if (use.latlong) {
        xpos <- sites[site, "long"]
        ypos <- sites[site, "lat"]
        do.call("floating.pie", c(list(xpos=xpos, ypos=ypos, x=xvals, radius=radius),
                            fill.args, ...))

      } else {
        ttl <- paste(sep="", site, " (n=", N, ")")
        do.call("pie", c(list(x=xvals, radius=radius, main=ttl, labels=""),
                            fill.args, ...))
      }
    }
    if (!is.na(add.legend) && !is.null(add.legend) && suppressWarnings(add.legend != "")) {
      ncol = floor(sqrt(length(fill.colors))) - 1
      n.per.code = apply(tab, 1, sum)
      legend.text = paste(sep="", names(fill.colors), " (", n.per.code, ")")
      legend(add.legend, ncol=ncol, legend=legend.text, col=fill.colors, 
             pch=15, bty="n", title="Code colors (N per code)",
             cex=1.0, x.intersp = 0.5, xpd=NA)
    }

  }
}


