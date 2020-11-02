whiten <- function (dens, lambda = 0.9) 
{
    dens$y <- lambda * dens$y + (1 - lambda)/diff(range(dens$x))
    dens
}

mysplit.density.plot <- function (obj, imp.vars = NULL, imp.vars.names = imp.vars, leg.posn = "topright", mfrow, mar = c(4.5, 1.5, 0.5, 4.5), omi = c(0.1, 0.25, 0.1, 0.1), bin = F, nbin = 101, leg.panel = 1, barwidth = 1, cex.legend = 0.8, line.ylab = 1.5, ...) 
{
    if (is.null(imp.vars)) 
        imp.vars <- imp.var.names <- names(sort(rowMeans(obj$imp.rsq), 
                                                decreasing = T))[1:2]
    is.binned <- function(obj) {
        compact <- obj$call$compact
        if (is.null(compact)) 
            FALSE
        else eval(compact)
    }
    normalize.histogram <- function(ci, integral = 1, bin = F, 
                                    nbin = 101) {
        if (bin) {
            brks <- seq(min(ci$x), max(ci$x), len = nbin)
            xx <- cut(ci$x, breaks = brks, inc = T)
            yy <- tapply(ci$y, xx, sum)
            yy[is.na(yy)] <- 0
            ci <- list(x = 0.5 * (brks[-1] + brks[-nbin]), y = yy)
        }
        dx <- min(diff(ci$x))
        Id <- sum(ci$y * dx)
        ci$y <- ci$y/Id * integral
        ci
    }
    normalize.density <- function(d, integral = 1, integrate = T) {
        Id <- if (integrate) 
            integrate.density(d)
        else 1
        d$y <- d$y/Id * integral
        d
    }
    integrate.density <- function(d) {
        integrate(approxfun(d, rule = 2), lower = min(d$x), upper = max(d$x))$value
    }
    scale.density <- function(d, scale = 1/mean(d$y)) {
        d$y <- d$y * scale
        d
    }
    par(mfrow = mfrow)
    nice.names <- structure(as.list(imp.vars.names), names = imp.vars)
    to_return <- vector(mode = "list", length = length(imp.vars))
    names(to_return) = imp.vars
    for (i in imp.vars) {
        imp <- importance(obj)[i]
        resA <- obj$res[obj$res$var == i, ]
        splits <- resA$split
        w <- pmax(resA$improve.norm, 0)
        X <- na.omit(obj$X[, i])
        rX <- range(X)
        dX <- diff(rX)
        dImp <- density(splits, weight = w/sum(w), from = rX[1], 
                        to = rX[2])
        if ((dX/dImp$bw) > 50) 
            dImp <- density(splits, weight = w/sum(w), from = rX[1], 
                            to = rX[2], bw = dX/50)
        dImpNorm <- normalize.density(dImp, imp, integrate = T)
        dObs <- density(X, from = rX[1], to = rX[2])
        if ((dX/dObs$bw) > 50) 
            dObs <- density(X, from = rX[1], to = rX[2], bw = dX/50)
        dObs <- whiten(dObs, lambda = 0.9)
        dObsNorm <- normalize.density(dObs, imp, integrate = T)
        ci <- cumimp(obj, i, standardize = F)
        ci$y <- diff(c(0, ci$y))
        ci <- normalize.histogram(ci, imp, bin = bin | !is.binned(obj), 
                                  nbin = nbin)
        dStd <- dImp
        dStd$y <- dImp$y/dObs$y
        dStdNorm <- try(normalize.density(dStd, imp, integrate = T))
        if (class(dStdNorm) == "try-error") 
            dStdNorm <- normalize.histogram(dStd, imp)
        plot(ci, type = "h", col = "grey60", xlim = range(splits), 
             lwd = barwidth, ylim = c(0, max(dImpNorm$y, dObsNorm$y, 
                                             dStdNorm$y) * 1.1), lend = 2, xlab = nice.names[[i]], 
             ylab = "", ...)
        lines(dImpNorm, col = "black", lwd = 2)
        lines(dObsNorm, col = "red", lwd = 2)
        lines(dStdNorm, col = "blue", lwd = 2)
        abline(h = mean(dStdNorm$y)/mean(dStd$y), lty = 2, col = "blue")
        if (i == imp.vars[leg.panel]) 
            legend(leg.posn, legend = c("Density of splits", 
                                        "Density of data", "Ratio of densities", "Ratio=1"), 
                   lty = c(1, 1, 1, 2), col = c("black", "red", 
                                                "blue", "blue"), cex = cex.legend, bty = "n", 
                   lwd = 1)
        myline = mean(dStdNorm$y)/mean(dStd$y)
        myreturn = list(dStdNorm, myline)
        names(myreturn) <- c("dStdNorm", "myline")
        to_return[[i]] <- myreturn
        }
    mtext("Density", side = 2, line = line.ylab, outer = T)
    return(to_return)
}

# given a returned density object from above, find local maxima and minima 
local_maxi <- function(mystdnorm) {
    require(pastecs)
    #make it a time series
    ts_y<-ts(mystdnorm$dStdNorm$y)
    #calculate turning points (extrema)
    tp<-turnpoints(ts_y)
    # get peaks 
    tpdf <- cbind(mystdnorm$dStdNorm$x[tp$peaks],
                  mystdnorm$dStdNorm$y[tp$peaks],
                  mystdnorm$dStdNorm$y[tp$peaks] > mystdnorm$myline
                  ) %>% 
        as.data.frame()
    colnames(tpdf) <- c("x", "y", "Acut")
    # #plot
    # plot(mystdnorm$dStdNorm)
    # points(tpdf[tpdf$Acut == T,'x'],
    #        tpdf[tpdf$Acut == T,'y'],
    #        col="red")
    return(tpdf)
}
