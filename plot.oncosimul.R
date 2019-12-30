A single object matching ‘plot.oncosimul’ was found
It was found in the following places
  registered S3 method for plot from namespace OncoSimulR
  namespace:OncoSimulR
with value

function (x, show = "drivers", type = ifelse(show == "genotypes", 
    "stacked", "line"), col = "auto", log = ifelse(type == "line", 
    "y", ""), ltyClone = 2:6, lwdClone = 0.9, ltyDrivers = 1, 
    lwdDrivers = 3, xlab = "Time units", ylab = "Number of cells", 
    plotClones = TRUE, plotDrivers = TRUE, addtot = FALSE, addtotlwd = 0.5, 
    ylim = NULL, xlim = NULL, thinData = FALSE, thinData.keep = 0.1, 
    thinData.min = 2, plotDiversity = FALSE, order.method = "as.is", 
    stream.center = TRUE, stream.frac.rand = 0.01, stream.spar = 0.2, 
    border = NULL, lwdStackedStream = 1, srange = c(0.4, 1), 
    vrange = c(0.8, 1), breakSortColors = "oe", legend.ncols = "auto", 
    ...) 
{
    if (!(type %in% c("stacked", "stream", "line"))) 
        stop("Type of plot unknown: it must be one of", "stacked, stream or line")
    if (!(show %in% c("genotypes", "drivers"))) 
        stop("show must be one of ", "genotypes or drivers")
    if (!(breakSortColors %in% c("oe", "distave", "random"))) 
        stop("breakSortColors must be one of ", "oe, distave, or random")
    colauto <- FALSE
    if ((length(col) == 1) && (col == "auto") && (type == "line") && 
        (show == "drivers")) 
        col <- c(8, "orange", 6:1)
    if ((length(col) == 1) && (col == "auto") && (show == "genotypes")) {
        col <- colorRampPalette(brewer.pal(8, "Dark2"))(ncol(x$pops.by.time) - 
            1)
        colauto <- TRUE
    }
    if (show == "genotypes") {
        plotDrivers <- FALSE
        plotClones <- TRUE
    }
    if (thinData) 
        x <- thin.pop.data(x, keep = thinData.keep, min.keep = thinData.min)
    if (!is.null(xlim)) 
        x <- xlim.pop.data(x, xlim)
    if (show == "drivers") {
        if (!inherits(x, "oncosimul2")) 
            ndr <- colSums(x$Genotypes[1:x$NumDrivers, , drop = FALSE])
        else {
            ndr <- colSums(x$Genotypes[x$Drivers, , drop = FALSE])
        }
    }
    else {
        ndr <- 1:(ncol(x$pops.by.time) - 1)
    }
    if ((type == "line") && is.null(ylim)) {
        if (log %in% c("y", "xy", "yx")) 
            ylim <- c(1, max(apply(x$pops.by.time[, -1, drop = FALSE], 
                1, sum)))
        else ylim <- c(0, max(apply(x$pops.by.time[, -1, drop = FALSE], 
            1, sum)))
    }
    if (plotDiversity) {
        oppd <- par(fig = c(0, 1, 0.8, 1))
        m1 <- par()$mar
        m <- m1
        m[c(1, 3)] <- c(0, 0.7)
        op <- par(mar = m)
        plotShannon(x)
        par(op)
        m1[c(3)] <- 0.2
        op <- par(mar = m1)
        par(fig = c(0, 1, 0, 0.8), new = TRUE)
    }
    if (plotClones) {
        plotClonesSt(x, ndr = ndr, show = show, na.subs = TRUE, 
            log = log, lwd = lwdClone, lty = ifelse(show == "drivers", 
                ltyClone, ltyDrivers), col = col, order.method = order.method, 
            stream.center = stream.center, stream.frac.rand = stream.frac.rand, 
            stream.spar = stream.spar, border = border, srange = srange, 
            vrange = vrange, type = type, breakSortColors = breakSortColors, 
            colauto = colauto, legend.ncols = legend.ncols, lwdStackedStream = lwdStackedStream, 
            xlab = xlab, ylab = ylab, ylim = ylim, xlim = xlim, 
            ...)
    }
    if (plotClones && plotDrivers && (type == "line")) 
        par(new = TRUE)
    if (plotDrivers && (type == "line")) {
        plotDrivers0(x, ndr, timescale = 1, trim.no.drivers = FALSE, 
            xlab = "", ylab = "", lwd = lwdDrivers, lty = ltyDrivers, 
            col = col, addtot = addtot, addtotlwd = addtotlwd, 
            log = log, ylim = ylim, xlim = xlim, legend.ncols = legend.ncols, 
            ...)
    }
    if (plotDiversity) {
        par(oppd)
    }
}
<bytecode: 0x5610a1da15f0>
<environment: namespace:OncoSimulR>
