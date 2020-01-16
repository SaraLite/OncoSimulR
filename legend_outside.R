##########################################################
############## Modified plot.oncosimul ###################
## The legend_outside.R script expands the possible parameters for the 
## representation of the fitness dependent frequency graphs. The new parameters 
## (legend.out and legend.pos) are only available when the script is loaded.

## legend.out -> Change the position of the legend. By default, the position 
## of the legend is on the top-left corner of the plot. If \code{TRUE},
## the legend is displaced outside the plot, on the rigth side. The new
## position is optimizied to Rmarkdown files. In the case that the canvas
## has a different size, see argument legend.out

## legend.pos -> If legend.out == TRUE, the position of the legend
## can be changed on the x axis by this argument. It allows to set a 
## inset value (argument of legend() function) with respect to 
## the rigth side of the plot.


## Function used by plot.oncosimul. Not modified
xlim.pop.data <- function(x, xlim) {
  x$pops.by.time <- x$pops.by.time[
    (x$pops.by.time[, 1] >= xlim[1]) &
      (x$pops.by.time[, 1] <= xlim[2]),   ]
  return(x)
}

## Function used by plot.oncosimul. Not modified
myhsvcols <- function(ndr, ymax, srange = c(0.4, 1),
                      vrange = c(0.8, 1),
                      breakSortColors = "oe") {
  ## Generate a set of colors so that:
  ##  - easy to tell when we increase number of drivers
  ##  - reasonably easy to tell in a legend
  ##  - different clones with same number of drivers have "similar" colors
  
  ## I use hsv color specification as this seems the most reasonable.
  
  minor <- table(ndr)
  major <- length(unique(ndr)) ## yeah same as length(minor), but least
  ## surprise
  
  h <- seq(from = 0, to = 1, length.out = major + 1)[-1]
  ## do not keep similar hues next to each other
  if(breakSortColors == "oe") {
    oe <- seq_along(h) %% 2
    h <- h[order(oe, h)]
  } else if(breakSortColors == "distave"){
    sl <- seq_along(h)
    h <- h[order(-abs(mean(sl) - sl))]
  } else if(breakSortColors == "random") {
    rr <- order(runif(length(h)))
    h <- h[rr]
  } 
  
  hh <- rep(h, minor)
  
  sr <- unlist(lapply(minor, function(x) 
    seq(from = srange[1], to = srange[2], length.out = x)))
  sv <- unlist(lapply(minor, function(x) 
    seq(from = vrange[1], to = vrange[2], length.out = x))
  )
  
  colors <- hsv(hh, sr, sv)
  
  ## This gives "average" or "median" color for legend
  ## colorsLegend <- aggregate(list(Color = colors), list(Drivers = ndr),
  ##                           function(x)
  ##                               as.character(x[((length(x) %/% 2) + 1 )]))
  
  ## Give the most abundant class color as the legend. Simpler to read
  colorsLegend <- by(data.frame(Color = colors, maxnum = ymax),
                     list(Drivers = ndr),
                     function(x) as.character(x$Color[which.max(x$maxnum)]))
  colorsLegend <- data.frame(Drivers = as.integer(row.names(colorsLegend)),
                             Color = cbind(colorsLegend)[, 1],
                             stringsAsFactors = FALSE)
  ## To show what it would look like
  ## plot(1:(sum(minor)), col = colors, pch = 16, cex = 3)
  ## legend(1, length(ndr), col = colorsLegend$Color, legend = names(minor),
  ##        pch = 16)
  
  return(list(colors = colors,
              colorsLegend = colorsLegend))
}

## Function used by plot.oncosimul. Not modified
plot.stacked2 <- function (x, y, order.method = "as.is", ylab = "", 
                           xlab = "", border = NULL, lwd = 1,
                           col = rainbow(length(y[1, ])), ylim = NULL, log = "", ...) 
{
  if (sum(y < 0) > 0) 
    stop("y cannot contain negative numbers")
  if (is.null(border)) 
    border <- par("fg")
  border <- as.vector(matrix(border, nrow = ncol(y), ncol = 1))
  col <- as.vector(matrix(col, nrow = ncol(y), ncol = 1))
  lwd <- as.vector(matrix(lwd, nrow = ncol(y), ncol = 1))
  if (order.method == "max") {
    ord <- order(apply(y, 2, which.max))
    y <- y[, ord]
    col <- col[ord]
    border <- border[ord]
  }
  if (order.method == "first") {
    ord <- order(apply(y, 2, function(x) min(which(x > 0))))
    y <- y[, ord]
    col <- col[ord]
    border <- border[ord]
  }
  top.old <- x * 0
  polys <- vector(mode = "list", ncol(y))
  for (i in seq(polys)) {
    top.new <- top.old + y[, i]
    polys[[i]] <- list(x = c(x, rev(x)), y = c(top.old, rev(top.new)))
    top.old <- top.new
  }
  if (is.null(ylim)) 
    ylim <- range(sapply(polys, function(x) range(x$y, na.rm = TRUE)), 
                  na.rm = TRUE)
  if (grepl("x", log)) 
    axes <- FALSE
  else axes <- TRUE
  plot(x, y[, 1], ylab = ylab, xlab = xlab, ylim = ylim, t = "n", 
       axes = axes, ...)
  for (i in seq(polys)) {
    polygon(polys[[i]], border = border[i], col = col[i], 
            lwd = lwd[i])
  }
  if (!axes) {
    relabelLogaxis(1)
    axis(2)
  }
}

## Main function modified at "line" and "stacked" legends plot
plot.oncosimul <- function(x,
                           show = "drivers", 
                           type = ifelse(show == "genotypes",
                                         "stacked", "line"),
                           col = "auto",
                           log = ifelse(type == "line", "y", ""),
                           ltyClone = 2:6,
                           lwdClone = 0.9,
                           ltyDrivers = 1,
                           lwdDrivers = 3,
                           xlab = "Time units",
                           ylab = "Number of cells",
                           plotClones = TRUE,
                           plotDrivers = TRUE,
                           addtot = FALSE,
                           addtotlwd = 0.5,
                           ylim = NULL,
                           xlim = NULL,
                           thinData = FALSE,
                           thinData.keep = 0.1,
                           thinData.min = 2,
                           plotDiversity = FALSE,
                           order.method = "as.is",
                           stream.center = TRUE,
                           stream.frac.rand = 0.01,
                           stream.spar = 0.2,
                           border = NULL,
                           lwdStackedStream = 1,
                           srange = c(0.4, 1),
                           vrange = c(0.8, 1),
                           breakSortColors = "oe",
                           legend.ncols = "auto",
                           legend.out = FALSE,
                           legend.pos = 0,
                           ...
) {
  
  
  if(!(type %in% c("stacked", "stream", "line")))
    stop("Type of plot unknown: it must be one of",
         "stacked, stream or line")
  
  if(!(show %in% c("genotypes", "drivers")))
    stop("show must be one of ",
         "genotypes or drivers")
  
  if(!(breakSortColors %in% c("oe", "distave", "random")))
    stop("breakSortColors must be one of ",
         "oe, distave, or random")
  
  
  
  colauto <- FALSE
  if((length(col) ==1) && (col == "auto") &&
     (type == "line") && (show == "drivers"))
    col <- c(8, "orange", 6:1)
  if((length(col) ==1) && (col == "auto") &&
     (show == "genotypes")) {
    ## For categorical data, I find Dark2, Paired, or Set1 to work best.
    col <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(ncol(x$pops.by.time) - 1)
    colauto <- TRUE
  }
  
  if(show == "genotypes") {
    plotDrivers <- FALSE
    plotClones <- TRUE
  }
  
  if(thinData)
    x <- thin.pop.data(x, keep = thinData.keep, min.keep = thinData.min)
  
  if(!is.null(xlim))
    x <- xlim.pop.data(x, xlim)
  
  ## For genotypes, ndr is now the genotypes.  Actually, ndr is now just
  ## a sequence 1:(ncol(y) - 1)
  
  ## The user will want to change the colors, like a colorRamp, etc. Or
  ## rainbow.
  
  ## genotypes and line, always call plotDrivers0
  if(show == "drivers") {
    if(!inherits(x, "oncosimul2"))
      ndr <- colSums(x$Genotypes[1:x$NumDrivers, , drop = FALSE])
    else {
      ndr <- colSums(x$Genotypes[x$Drivers, , drop = FALSE])
    }
  } else { ## show we are showing genotypes
    ndr <- 1:(ncol(x$pops.by.time) - 1)
  }
  
  if((type == "line") && is.null(ylim)) {
    if(log %in% c("y", "xy", "yx") )
      ylim <- c(1, max(apply(x$pops.by.time[, -1, drop = FALSE], 1, sum)))
    else
      ylim <- c(0, max(apply(x$pops.by.time[, -1, drop = FALSE], 1, sum)))
  }
  if(plotDiversity) {
    oppd <- par(fig = c(0, 1, 0.8, 1))
    m1 <- par()$mar
    m <- m1
    m[c(1, 3)] <- c(0, 0.7)
    op <- par(mar = m )
    plotShannon(x)
    par(op)
    m1[c(3)] <- 0.2
    op <- par(mar = m1)
    par(fig = c(0, 1, 0, 0.8), new = TRUE)  
  }
  
  ## Shows its history: plotClones makes plotDrivers0 unneeded with
  ## stacked and stream plots. But now so with line plot.
  ## When showing genotypes, plotDrivers0 with line only used for
  ## showing the legend.
  if(plotClones) {
    plotClonesSt(x,
                 ndr = ndr,
                 show = show,
                 na.subs = TRUE,
                 log = log,
                 lwd = lwdClone,
                 lty = ifelse(show == "drivers", ltyClone, ltyDrivers),
                 col = col, 
                 order.method = order.method,
                 stream.center = stream.center,
                 stream.frac.rand = stream.frac.rand,
                 stream.spar = stream.spar,
                 border = border,
                 srange = srange,
                 vrange = vrange,
                 type = type,
                 breakSortColors = breakSortColors,
                 colauto = colauto,
                 legend.ncols = legend.ncols,
                 lwdStackedStream = lwdStackedStream,
                 xlab = xlab,
                 ylab = ylab,
                 ylim = ylim,
                 xlim = xlim,
                 legend.out = legend.out,
                 legend.pos = legend.pos,
                 ...)
  }
  
  if(plotClones && plotDrivers && (type == "line"))
    par(new = TRUE)
  
  if( plotDrivers && (type == "line") ) {
    plotDrivers0(x,
                 ndr,
                 timescale = 1,
                 trim.no.drivers = FALSE,
                 xlab = "", ylab = "",
                 lwd = lwdDrivers,
                 lty = ltyDrivers,
                 col = col, 
                 addtot = addtot,
                 addtotlwd = addtotlwd,
                 log = log, ylim = ylim,
                 xlim = xlim,
                 legend.ncols = legend.ncols,
                 ...)
  }
  if(plotDiversity) {
    par(oppd)
  }
  
}

plotClonesSt <- function(z,
                         ndr,
                         show = "drivers",
                         na.subs = TRUE,
                         log = "y",
                         lwd = 1,
                         ## type = "l",
                         lty = 1:8, col = 1:9,
                         order.method = "as.is",
                         stream.center = TRUE,
                         stream.frac.rand = 0.01,
                         stream.spar = 0.2,
                         border = NULL,
                         srange = c(0.4, 1),
                         vrange = c(0.8, 1),
                         type = "stacked",
                         breakSortColors = "oe",
                         colauto = TRUE,
                         legend.ncols = "auto",
                         lwdStackedStream = 1,
                         xlab = "Time units",
                         ylab = "Number of cells",
                         ylim = NULL,
                         xlim = NULL,
                         legend.out = FALSE,
                         legend.pos = 0,
                         ...) {
  
  ## if given ndr, we order columns based on ndr, so clones with more
  ## drivers are plotted last
  
  y <- z$pops.by.time[, 2:ncol(z$pops.by.time), drop = FALSE]
  
  ## Code in stacked and stream plots relies on there being no NAs. Could
  ## change it, but it does not seem reasonable.
  ##  But my original plotting code runs faster and is simpler if 0 are
  ##  dealt as NAs (which also makes log transformations simpler).
  
  if(type %in% c("stacked", "stream") )
    na.subs <- FALSE
  
  if(na.subs){
    y[y == 0] <- NA
  }
  ## if(is.null(ndr))
  ##     stop("Should never have null ndr")
  ## if(!is.null(ndr)) {
  ## could be done above, to avoid creating
  ## more copies
  oo <- order(ndr)
  y <- y[, oo, drop = FALSE]
  ndr <- ndr[oo]
  if(show == "drivers") {
    col <- rep(col, length.out = (1 + max(ndr)))[ndr + 1]
    lty <- rep(lty, length.out = ncol(y))
  } else {
    if(length(col) < max(ndr))
      warning("Repeating colors; you might want to",
              "pass a col vector of more elements")
    col <- rep(col, length.out = (max(ndr)))[ndr]
  }
  
  ## Set variables to place the legend if legend.out = TRUE
  if (legend.out) {
    mar = c(4, 4.8, 3, 6)
    xpd = TRUE
    position = "right"
    ## If legend.pos is zero (by default), the inset value is -0.2
    ## This value allows a correct position in Rmarkdown files
    if (legend.pos == 0)
      inset = -0.2
    else
      inset = legend.pos
  } else {
    mar = c(5, 4, 4, 2) + 0.1
    xpd = FALSE
    position = "topleft"
    inset = 0
  }
  
  if(type == "line") {
    par(mar = mar)
    matplot(x = z$pops.by.time[, 1],
            y = y,
            log = log, type = "l",
            col = col, lty = lty,
            lwd = lwd,
            xlab = xlab,
            ylab = ylab,
            ylim = ylim,
            xlim = xlim,
            ...)
    box()
    if(show == "genotypes") {
      if(!inherits(z, "oncosimul2")) {
        ldrv <- genotypeLabel(z)
      } else {
        ldrv <- z$GenotypesLabels
      }
      ldrv[ldrv == ""] <- "WT"
      ldrv[ldrv == " _ "] <- "WT"
      if(legend.ncols == "auto") {
        if(length(ldrv) > 6) legend.ncols <- 2
        else legend.ncols <- 1
      }
      par(xpd = xpd)
      legend(x = position, title = "Genotypes", lty = lty, 
             inset = inset, col = col, lwd = lwd, legend = ldrv, 
             ncol = legend.ncols)
    }
  } else {
    par(mar = mar)
    ymax <- colSums(y)
    if((show == "drivers") || ((show == "genotypes") && (colauto))) {
      cll <- myhsvcols(ndr, ymax, srange = srange, vrange = vrange,
                       breakSortColors = breakSortColors)
    } else {
      cll <- list(colors = col)
    }
    x <- z$pops.by.time[, 1]
    if(grepl("y", log)) {
      stop("It makes little sense to do a stacked/stream",
           "plot after taking the log of the y data.")
    }
    if(grepl("x", log)) {
      x <- log10(x + 1)
    }
    
    if (type == "stacked") {
      plot.stacked2(x = x,
                    y = y,
                    order.method = order.method,
                    border = border,
                    lwd = lwdStackedStream,
                    col = cll$colors,
                    log = log,
                    xlab = xlab,
                    ylab = ylab,
                    ylim = ylim,
                    xlim = xlim,
                    ...) 
    } else if (type == "stream") {
      plot.stream2(x = x,
                   y = y,
                   order.method = order.method,
                   border = border,
                   lwd = lwdStackedStream,
                   col = cll$colors,
                   frac.rand = stream.frac.rand,
                   spar = stream.spar,
                   center = stream.center,
                   log = log,
                   xlab = xlab,
                   ylab = ylab,
                   ylim = ylim,
                   xlim = xlim,
                   ...)
    }
    if(show == "drivers") {
      par(mar = mar)
      if(legend.ncols == "auto") {
        if(length(cll$colorsLegend$Drivers) > 6) legend.ncols <- 2
        else legend.ncols <- 1
      }
      par(xpd = xpd)
      legend(x = position, title = "Number of drivers", 
             inset = inset, pch = 15, 
             col = cll$colorsLegend$Color, 
             legend = cll$colorsLegend$Drivers, ncol = legend.ncols)
      
    } else if (show == "genotypes") {
      par(mar = mar)
      if(!inherits(z, "oncosimul2")) {
        ldrv <- genotypeLabel(z)
      } else {
        ldrv <- z$GenotypesLabels
      }
      ldrv[ldrv == ""] <- "WT"
      ldrv[ldrv == " _ "] <- "WT"            
      if(legend.ncols == "auto") {
        if(length(ldrv) > 6) legend.ncols <- 2
        else legend.ncols <- 1
      }
      par(xpd = xpd)
      legend(x = position, title = "Genotypes", pch = 15,
             inset = inset, lty = lty, lwd = lwd,
             col = cll$colors, legend = ldrv, ncol = legend.ncols)
    }
  }
}
