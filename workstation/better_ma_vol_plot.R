## Draw MA / Volcano plot
drawMA <- function(resultDF, qco, fco, ttl_pre = "", ylim = NULL, xlim = NULL, mode = "ma", pcol = "padj", ...) {

	if(!mode %in% c("ma", "vol")) stop("Argument mode should be either 'ma' or 'vol'.")
	if(!pcol %in% c("pval", "padj")) stop("Argument pcol should be either 'pval' or 'padj'.")
	if(pcol == "pval") { mt_p <- "pv"; ax_p <- "-log10(pval)" }
	if(pcol == "padj") { mt_p <- "qv"; ax_p <- "-log10(padj)" }

	resultDF <- resultDF %>% mutate(p.used = get(pcol)) 

	ui <- which(resultDF$p.used < qco & resultDF$logFC >=  log2(fco))
	di <- which(resultDF$p.used < qco & resultDF$logFC <= -log2(fco))
	gcnt <- paste("up:", length(ui), "dn:", length(di))
	mtitle <- paste(ttl_pre, "fc", fco, mt_p, qco, gcnt)
	lfcmx <- max(abs(resultDF$logFC), na.rm = TRUE)
	axlim <- c(-lfcmx, lfcmx)

	# MA
	if( mode == "ma" ) {
		plot(logFC ~ AveExpr, data = resultDF, pch = 20, main = mtitle, cex = 0.05, ylim = axlim, ...)
		points(logFC ~ AveExpr, data = resultDF[c(ui, di), ], pch = 20, col = "red", cex = 0.25)
		# abline(h = 0, col = "blue", lty = 2)
		abline(h = c(-log2(fco), log2(fco)), col = "blue", lty = 2)
	}
	
	# volcano
	if( mode == "vol" ) {
		plot(-log10(p.used) ~ logFC, data = resultDF, pch = 20, main = mtitle, cex = 0.05, xlim = axlim, ylab = ax_p, ...)
		points(-log10(p.used) ~ logFC, data = resultDF[c(ui, di), ], pch = 20, col = "red", cex = 0.25)
		# abline(v = 0, col = "blue", lty = 2)
		abline(v = c(-log2(fco), log2(fco)), col = "blue", lty = 2)
		abline(h = -log10(qco), col = "blue", lty = 2)
	}
}



drawMA2 <- function(resultDF, qco, fco, ttl_pre = "", ylim = NULL, xlim = NULL, mode = "ma", pcol = "padj", top10_lab = FALSE) {
	require(ggrepel)

	if(!mode %in% c("ma", "vol")   ) stop("Argument mode should be either 'ma' or 'vol'."   )
	if(!pcol %in% c("pval", "padj")) stop("Argument pcol should be either 'pval' or 'padj'.")
	if(pcol == "pval") { mt_p <- "pv"; ax_p <- expression(bold(-log[10] ~ p.val)) }
	if(pcol == "padj") { mt_p <- "qv"; ax_p <- expression(bold(-log[10] ~ FDR))   }

	resultDF <- resultDF %>% mutate(p.used = get(pcol)) %>% filter(!is.na(p.used))

	ui <- which(resultDF$p.used < qco & resultDF$logFC >=  log2(fco))
	di <- which(resultDF$p.used < qco & resultDF$logFC <= -log2(fco))
	gcnt <- paste("up:", length(ui), "dn:", length(di))
	mtitle <- paste(ttl_pre, "FC", fco, mt_p, qco, gcnt)
	lfcmx <- max(abs(resultDF$logFC), na.rm = TRUE)
	axlim <- c(-lfcmx, lfcmx)

	# if(diff(round(axlim)) >= 10) {
	# 	x_ax <- seq(round(axlim/2)[1], round(axlim/2)[2], 1) * 2
	# } else {
	# 	x_ax <- seq(round(axlim)[1], round(axlim)[2], 1)
	# }

	# MA
	if( mode == "ma" ) {
		gp <- ggplot(resultDF, aes(x = AveExpr, y = logFC)) +
			geom_point(shape = 20, cex = 0.05) +
			geom_point(data = resultDF[c(ui, di), ], aes(x = AveExpr, y = logFC), shape = 20, colour = "red", cex = 0.25) +
			geom_hline(yintercept = c(-log2(fco), log2(fco)), col = "blue", lty = 2) +
			scale_y_continuous(limits = axlim) + 
			# scale_y_continuous(breaks = x_ax, limits = axlim) + 
			labs(title = mtitle, y = expression(bold(log[2] ~ FC)), x = "AveExpr")
	}
	
	# volcano
	if( mode == "vol" ) {
		gp <- ggplot(resultDF, aes(x = logFC, y = -log10(p.used))) +
			geom_point(shape = 20, cex = 0.05) +
			geom_point(data = resultDF[c(ui, di), ], aes(x = logFC, y = -log10(p.used)), shape = 20, colour = "red", cex = 0.25) +
			geom_hline(yintercept = -log10(qco), col = "blue", lty = 2) + 
			geom_vline(xintercept = c(-log2(fco), log2(fco)), col = "blue", lty = 2) +
			scale_x_continuous(limits = axlim) +
			# scale_x_continuous(breaks = x_ax, limits = axlim) +
			labs(title = mtitle, y = ax_p, x = expression(bold(log[2] ~ FC)))
	}

	# label
	if( top10_lab ) {
		labeldf <- resultDF[c(ui, di), ] %>% 
			filter(rank(logFC, ties.method = "min") <= 10 | rank(-logFC, ties.method = "min") <= 10)
		gp <- gp + 
			geom_label_repel(data = labeldf, aes(x = logFC, y = -log10(p.used), label = geneSym), colour = "blue", max.overlaps = 20)
	}

	gp <- gp + 
		theme(
			panel.border = element_blank(),
			panel.grid = element_line(colour = "grey92"), 
			plot.title = element_text(colour = "black", face = "bold", size = rel(1.4), hjust = 0.5),
			axis.title = element_text(colour = "black", face = "bold", size = rel(1.1)),
			axis.line  = element_line(colour = "black", linewidth = 0.5),
			axis.text  = element_text(colour = "black", size = rel(1.1)),
			plot.background   = element_rect(fill = "transparent", colour = NA), 
			panel.background  = element_rect(fill = "transparent", colour = NA), 
			strip.background  = element_rect(fill = "transparent", colour = NA, linewidth = 0.7), 
			legend.background = element_rect(fill = "transparent", colour = NA), 
			legend.key = element_blank(), 
			panel.ontop = FALSE,
			complete = TRUE
		)

	gp
}
