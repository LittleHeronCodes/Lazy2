
##====================================================
## GSEA plot

gset = intgpath[[i]]
stats = tstatLS$PHMG
nes=ss1$NES
qv=ss1$qVal
mtitle='PHMG'
ylab=cleanTermNames(i, remove_source=TRUE)
lwd=1.1
base_size=11
ticksSize=0.3
ylim = c(-0.55,0.55)
line.col='green'
draw=TRUE

plotEnrichment2 <- function(gset, stats, nes, qv, gseaParam = 1, mtitle=NULL, ylab='',
	ticksSize=0.4, base_size=7, line.col='green', lwd=2, ylim=NULL, draw=TRUE, statbar=FALSE) {

	require(grid)
	require(gtable)
	require(fgsea)

	rnk <- rank(-stats)
	ord <- order(rnk)
	stats <- stats[ord]
	statsAdj <- sign(stats) * (abs(stats)^gseaParam)
	statsAdj <- statsAdj/max(abs(statsAdj))
	pathway <- unname(as.vector(na.omit(match(gset, names(statsAdj)))))
	pathway <- sort(pathway)
	gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, returnAllExtremes = TRUE)
	bottoms <- gseaRes$bottoms
	tops <- gseaRes$tops

	## Add fgsea run if no nes and qv given ##
	if(any(is.na(c(nes,qv)))) {
		
	}

	txt <- sprintf('NES : %.2f  \nq-value : %.2e  ', nes, qv)

	n <- length(statsAdj)
	xs <- as.vector(rbind(pathway - 1, pathway))
	ys <- as.vector(rbind(bottoms, tops))
	toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
	diff <- (max(tops) - min(bottoms))/8
	es <- c(tops, bottoms)[which.max(abs(c(tops,bottoms)))]

	if(is.null(ylim)) ylim = c(floor(min(bottoms)*20), ceiling(max(tops)*20))/20

	if( (ylim[2] - ylim[1]) > 0.55 ) {
		br <- 0.2
		ln1 <- seq(ceiling(ylim[1]*5)/5, floor(ylim[2]*5)/5, br)
	} else {
		br <- 0.1
		ln1 <- seq(ceiling(ylim[1]*10)/10, floor(ylim[2]*10)/10, br)
	}

	txt_y_pos = ifelse(es < 0, max(c(ln1,tops)), min(c(ln1, bottoms)))

	half_line = base_size/2

	g1 <- ggplot(toPlot, aes(x = x, y = y)) + geom_point(color = line.col, size = 0.1) + 
	  geom_hline(yintercept = 0, colour = "black", linetype='dashed', size = lwd*0.8) +
	  geom_hline(yintercept = ifelse(nes>0, max(tops), min(bottoms)), colour = "red", linetype = "dashed", size = lwd*0.8) +
	  geom_line(color = line.col,size = lwd) +
	  annotate('text', x=max(toPlot$x), y=txt_y_pos, label=txt, hjust=1, vjust=.5, fontface='plain', size=rel(3.0)) +
	  scale_y_continuous(breaks=ln1, limits=ylim) + 
	  labs(y = ylab, title=mtitle) + 
	  theme_common_gsea(base_size=base_size) +
	  theme(
	  	panel.grid.major.y = element_line(colour = "grey85",linetype='dashed', size = lwd*0.65),
	    plot.title = element_text(hjust = 0.5, vjust=0.2, face='bold', margin=unit(c(0,0,1.5,0), 'mm')),
	    axis.title.y = element_text(face='bold', angle=90, margin=unit(c(0,1.5,0,0), 'mm'),size=rel(0.95)),
	    # plot.margin = unit(c(1.2,2.0,0,1.2), 'mm'),
	    plot.margin = margin(half_line, half_line*2.5, 0, half_line),
	  	axis.title.x=element_blank(), axis.text.x=element_blank()
	  	)

	g2 <- ggplot(data.frame(x=pathway),aes(x = x, y = -diff/3, xend = x, yend = diff/3)) +
	  geom_segment(size = ticksSize, colour='grey35') +
	  theme_common_gsea(base_size=base_size) + scale_x_continuous(breaks=c(0,5000,10000,15000), limits=c(0,15000)) +
	  theme(
		# plot.margin = unit(c(0,2.0,1.2,1.2), 'mm'),
		plot.margin = margin(0, half_line*2.5, half_line, half_line),
	  	axis.text.y = element_blank(), axis.title=element_blank()
	 	)

	if(statbar) {
		lead <- toPlot[which(toPlot$y==es),]
		g1 <- g1 + geom_segment(x=lead$x, xend=lead$x, y=0, yend=es, color='red', linetype=3, size=lwd*0.25)
		g2 <- g2 + theme(axis.text.x=element_blank(),plot.margin = margin(0, half_line*2.5, 0, half_line))
		g3 <- ggplot(data.frame(x=seq_along(stats), stat=stats), aes(x=x, y=0, xend=x, yend=stat)) +
		  geom_segment(size = ticksSize, colour='grey35') +
		  theme_common_gsea(base_size=base_size) +
		  theme(
		  	plot.margin = margin(0, half_line*2.5, half_line, half_line),
		  	axis.text.y=element_blank(), axis.title=element_blank()
		  	)
	}

	gr1 <- ggplotGrob(g1)
	gr2 <- ggplotGrob(g2)

	if(statbar) {
		gr3 <- ggplotGrob(g3)
		ht.ratio <- c(6,1,1)
		gr <- rbind(gr1, gr2, gr3)
	} else {
		ht.ratio <- c(7,1)
		gr <- rbind(gr1, gr2)		
	}

	# gr$widths <- grid::unit.pmax(gr1$widths, gr2$widths, gr3$widths)
	gr$widths <- grid::unit.pmax(gr1$widths, gr2$widths)

	# identify the position of the panels within the gtable
	panid <- gr$layout$t[grep(pattern="panel", gr$layout$name)]
	gr$heights[panid] <- unit(ht.ratio, 'null')

	if(draw) {
		grid.newpage()
		grid.draw(gr)		
	}

	return(gr)
}


theme_common_gsea <- function(base_size=5) {
	half_line = base_size/2
	.theme <- theme(
		text = element_text(face='plain', size=base_size, colour='black', family='Arial'),
		plot.title = element_text(size=rel(1.0)),
		axis.ticks=element_blank(), 
		axis.text = element_text(size=rel(0.8), colour='black'),
		panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		panel.background = element_rect(fill = "transparent", colour = NA),
		plot.background = element_rect(fill = "transparent",colour = NA),
		complete=TRUE)
	.theme
}

