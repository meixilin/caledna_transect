# change the species.cumulative.plot

myspecies.cumulative.plot <- function (obj, imp.vars = NULL, imp.vars.names = imp.vars, leg.posn = "topleft", leg.nspecies = 10, legend = TRUE, 
                                       mfrow, show.species = TRUE, 
                                       show.overall = TRUE, mar = c(0, 2.1, 1.1, 0), omi = c(0.75, 0.75, 0.1, 0.1), common.scale = F, 
                                       line.ylab = 1, cex.legend = 0.75, ...) 
{
    # modified by: Meixi Lin, source code: ~/UCLA/abiotic_factors/r_codes/step4_gradient_forest/change_species_cumulative_plot.R
    if (is.null(imp.vars)) 
        imp.vars <- imp.var.names <- names(importance(obj))[1:2]
    par(mfrow = mfrow)
    cols <- rainbow(length(levels(obj$res.u$spec)))
    names(cols) <- levels(obj$res.u$spec)
    xaxt <- if (show.overall) 
        "n"
    else "s"
    if (show.species) {
        for (varX in imp.vars) {
            CU <- cumimp(obj, varX, "Species")
            xlim <- range(sapply(CU, "[[", "x"), na.rm = T)
            ylim <- range(sapply(CU, "[[", "y"), na.rm = T)
            plot(xlim, ylim, type = "n", xlab = if (show.overall) 
                ""
                else imp.vars.names[imp.vars == varX], ylab = "", 
                xaxt = xaxt, ...)
            for (species in names(CU)) {
                isub <- seq(1, length(CU[[species]]$x), len = pmin(500, 
                                                                   length(CU[[species]]$x)))
                lines(CU[[species]]$x[isub], CU[[species]]$y[isub], 
                      type = "s", col = cols[species])
            }
            no.species <- length(names(cols))
            imp.sp <- sapply(CU, function(cu) max(cu$y))
            best <- order(-imp.sp)[1:min(leg.nspecies, length(imp.sp))]
            if (legend) 
                legend(x = leg.posn, legend = names(cols)[best], 
                       pch = rep(1, no.species)[best], col = cols[best], 
                       bty = "n", cex = cex.legend, pt.lwd = 2)
        }
    }
    if (show.overall) {
        for (varX in imp.vars) {
            CU <- cumimp(obj, varX)
            ymax <- max(CU$y, na.rm = T)
            if (varX == imp.vars[1]) 
                ymax1 <- ymax
            isub <- seq(1, length(CU$x), len = pmin(500, length(CU$x)))
            plot(CU$x[isub], CU$y[isub], type = "s", ylab = "", 
                 xlab = imp.vars.names[imp.vars == varX], ylim = c(0, 
                                                                   if (common.scale) ymax1 else ymax), ...)
        }
    }
    mtext("Cumulative importance", side = 2, line = line.ylab, 
          outer = TRUE)
}

# # Original archive --------
# function (obj, imp.vars = NULL, imp.vars.names = imp.vars, leg.posn = "topleft", 
#           leg.nspecies = 10, legend = TRUE, mfrow = rev(n2mfrow(length(imp.vars) * 
#                                                                     (show.species + show.overall))), show.species = TRUE, 
#           show.overall = TRUE, mar = c(0, 2.1, 1.1, 0), omi = c(0.75, 
#                                                                 0.75, 0.1, 0.1), common.scale = F, line.ylab = 1, cex.legend = 0.75, 
#           ...) 
# {
#     if (is.null(imp.vars)) 
#         imp.vars <- imp.var.names <- names(importance(obj))[1:2]
#     par(mfrow = mfrow)
#     cols <- rainbow(length(levels(obj$res.u$spec)))
#     names(cols) <- levels(obj$res.u$spec)
#     xaxt <- if (show.overall) 
#         "n"
#     else "s"
#     if (show.species) {
#         for (varX in imp.vars) {
#             CU <- cumimp(obj, varX, "Species")
#             xlim <- range(sapply(CU, "[[", "x"))
#             ylim <- range(sapply(CU, "[[", "y"))
#             plot(xlim, ylim, type = "n", xlab = if (show.overall) 
#                 ""
#                 else imp.vars.names[imp.vars == varX], ylab = "", 
#                 xaxt = xaxt, ...)
#             for (species in names(CU)) {
#                 isub <- seq(1, length(CU[[species]]$x), len = pmin(500, 
#                                                                    length(CU[[species]]$x)))
#                 lines(CU[[species]]$x[isub], CU[[species]]$y[isub], 
#                       type = "s", col = cols[species])
#             }
#             no.species <- length(names(cols))
#             imp.sp <- sapply(CU, function(cu) max(cu$y))
#             best <- order(-imp.sp)[1:min(leg.nspecies, length(imp.sp))]
#             if (legend) 
#                 legend(x = leg.posn, legend = names(cols)[best], 
#                        pch = rep(1, no.species)[best], col = cols[best], 
#                        bty = "n", cex = cex.legend, pt.lwd = 2)
#         }
#     }
#     if (show.overall) {
#         for (varX in imp.vars) {
#             CU <- cumimp(obj, varX)
#             ymax <- max(CU$y)
#             if (varX == imp.vars[1]) 
#                 ymax1 <- ymax
#             isub <- seq(1, length(CU$x), len = pmin(500, length(CU$x)))
#             plot(CU$x[isub], CU$y[isub], type = "s", ylab = "", 
#                  xlab = imp.vars.names[imp.vars == varX], ylim = c(0, 
#                                                                    if (common.scale) ymax1 else ymax), ...)
#         }
#     }
#     mtext("Cumulative importance", side = 2, line = line.ylab, 
#           outer = TRUE)
# }
# # <bytecode: 0x164a667f0>
# #     <environment: namespace:gradientForest>