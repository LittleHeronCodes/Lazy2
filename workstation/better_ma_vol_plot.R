
## Draw MA / Volcano plot one function
drawMA <- function(resultDF, qco, fco, ttl_pre="", ylim = NULL, xlim = NULL, mode = "ma") {
	if(!mode %in% c("ma", "vol")) stop("Argument mode should be either 'ma' or 'vol'.")

    ui <- which(resultDF$padj < qco & resultDF$logFC > log2(fco))
    di <- which(resultDF$padj < qco & resultDF$logFC < -log2(fco))
    gcnt <- paste("up:", length(ui), "dn:", length(di))

    if( mode == "ma" & is.null(ylim) ) {
        ylim <- c(-max(abs(resultDF$logFC), na.rm = TRUE), max(abs(resultDF$logFC), na.rm = TRUE))
    }
	if( mode == "vol" & is.null(xlim)) {
        xlim <- c(-max(abs(resultDF$logFC), na.rm = TRUE), max(abs(resultDF$logFC), na.rm = TRUE))
	}
    mtitle <- paste(ttl_pre, "fc", fco, "qv", qco, gcnt)

    # MA
	if( mode == "ma" ) {
		plot(logFC ~ AveExpr, data = resultDF, pch = 20, main = mtitle, cex = 0.05, ylim = ylim)
		points(logFC ~ AveExpr, data = resultDF[c(ui, di), ], pch = 20, col = "red", cex = 0.25)
		abline(h = 0, col = "blue", lty = 2)
	}
    # volcano
	if( mode == "vol" ) {
		plot(-log10(padj) ~ logFC, data = resultDF, pch = 20, main = mtitle, cex = 0.05, xlim = xlim)
		points(-log10(padj) ~ logFC, data = resultDF[c(ui, di), ], pch = 20, col = "red", cex = 0.25)
		abline(v = 0, col = "blue", lty = 2)
	}
}
