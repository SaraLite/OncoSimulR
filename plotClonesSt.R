A single object matching ‘plotClonesSt’ was found
It was found in the following places
  namespace:OncoSimulR
with value

function (z, ndr, show = "drivers", na.subs = TRUE, log = "y", 
    lwd = 1, lty = 1:8, col = 1:9, order.method = "as.is", stream.center = TRUE, 
    stream.frac.rand = 0.01, stream.spar = 0.2, border = NULL, 
    srange = c(0.4, 1), vrange = c(0.8, 1), type = "stacked", 
    breakSortColors = "oe", colauto = TRUE, legend.ncols = "auto", 
    lwdStackedStream = 1, xlab = "Time units", ylab = "Number of cells", 
    ylim = NULL, xlim = NULL, ...) 
{
    y <- z$pops.by.time[, 2:ncol(z$pops.by.time), drop = FALSE]
    if (type %in% c("stacked", "stream")) 
        na.subs <- FALSE
    if (na.subs) {
        y[y == 0] <- NA
    }
    oo <- order(ndr)
    y <- y[, oo, drop = FALSE]
    ndr <- ndr[oo]
    if (show == "drivers") {
        col <- rep(col, length.out = (1 + max(ndr)))[ndr + 1]
        lty <- rep(lty, length.out = ncol(y))
    }
    else {
        if (length(col) < max(ndr)) 
            warning("Repeating colors; you might want to", "pass a col vector of more elements")
        col <- rep(col, length.out = (max(ndr)))[ndr]
    }
    if (type == "line") {
        matplot(x = z$pops.by.time[, 1], y = y, log = log, type = "l", 
            col = col, lty = lty, lwd = lwd, xlab = xlab, ylab = ylab, 
            ylim = ylim, xlim = xlim, ...)
        box()
        if (show == "genotypes") {
            if (!inherits(z, "oncosimul2")) {
                ldrv <- genotypeLabel(z)
            }
            else {
                ldrv <- z$GenotypesLabels
            }
            ldrv[ldrv == ""] <- "WT"
            ldrv[ldrv == " _ "] <- "WT"
            if (legend.ncols == "auto") {
                if (length(ldrv) > 6) 
                  legend.ncols <- 2
                else legend.ncols <- 1
            }
            legend(x = "topleft", title = "Genotypes", lty = lty, 
                col = col, lwd = lwd, legend = ldrv, ncol = legend.ncols)
        }
    }
    else {
        ymax <- colSums(y)
        if ((show == "drivers") || ((show == "genotypes") && 
            (colauto))) {
            cll <- myhsvcols(ndr, ymax, srange = srange, vrange = vrange, 
                breakSortColors = breakSortColors)
        }
        else {
            cll <- list(colors = col)
        }
        x <- z$pops.by.time[, 1]
        if (grepl("y", log)) {
            stop("It makes little sense to do a stacked/stream", 
                "plot after taking the log of the y data.")
        }
        if (grepl("x", log)) {
            x <- log10(x + 1)
        }
        if (type == "stacked") {
            plot.stacked2(x = x, y = y, order.method = order.method, 
                border = border, lwd = lwdStackedStream, col = cll$colors, 
                log = log, xlab = xlab, ylab = ylab, ylim = ylim, 
                xlim = xlim, ...)
        }
        else if (type == "stream") {
            plot.stream2(x = x, y = y, order.method = order.method, 
                border = border, lwd = lwdStackedStream, col = cll$colors, 
                frac.rand = stream.frac.rand, spar = stream.spar, 
                center = stream.center, log = log, xlab = xlab, 
                ylab = ylab, ylim = ylim, xlim = xlim, ...)
        }
        if (show == "drivers") {
            if (legend.ncols == "auto") {
                if (length(cll$colorsLegend$Drivers) > 6) 
                  legend.ncols <- 2
                else legend.ncols <- 1
            }
            legend(x = "topleft", title = "Number of drivers", 
                pch = 15, col = cll$colorsLegend$Color, legend = cll$colorsLegend$Drivers, 
                ncol = legend.ncols)
        }
        else if (show == "genotypes") {
            if (!inherits(z, "oncosimul2")) {
                ldrv <- genotypeLabel(z)
            }
            else {
                ldrv <- z$GenotypesLabels
            }
            ldrv[ldrv == ""] <- "WT"
            ldrv[ldrv == " _ "] <- "WT"
            if (legend.ncols == "auto") {
                if (length(ldrv) > 6) 
                  legend.ncols <- 2
                else legend.ncols <- 1
            }
            legend(x = "topleft", title = "Genotypes", pch = 15, 
                col = cll$colors, legend = ldrv, ncol = legend.ncols)
        }
    }
}
<bytecode: 0x5610a8084d08>
<environment: namespace:OncoSimulR>
