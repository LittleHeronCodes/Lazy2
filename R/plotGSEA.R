#' plotEnrichment2
#'
#' Plot running score for GSEA enrichment plot. Function modified from plotEnrichment function in fgsea package.
#' @param gset gene set to calculate enrichment in character vector.
#' @param stats gene level statistics. This should be a named vector where names match gene ids in gset.
#' @param nes NES provided for annotation. Runs fgsea if not provided.
#' @param qv adjusted p value provided for annotation. Runs fgsea if not provided.
#' @param gseaParam gseaParam from plotEnrichment.
#' @param mtitle main title
#' @param ytitle y axis name
#' @param ticksSize ticks lwd for geneset bars.
#' @param base_size font base size. default 7
#' @param line.col line colour for GSEA running score.
#' @param lwd line width for GSEA running score.
#' @param ylims y-axis limits.
#' @param draw draw plot. default TRUE
#' @param statbar show statistics tick at bottom of plot.
#' @return enrichment plot
#' @import grid
#' @importFrom fgsea fgsea calcGseaStat
#' @export

## GSEA plot
plotEnrichment2 <- function(gset, stats, nes, qv, gseaParam = 1, mtitle = NULL, ytitle = "",
                            ticksSize = 0.4, base_size = 7, line.col = "green", lwd = 2, ylims = NULL, draw = TRUE, statbar = FALSE) {
    set.seed(1234)

    if (!any(gset %in% names(stats))) {
        stop("all genes in gset should be in names(stats).")
    }

    rnk <- rank(-stats)
    ord <- order(rnk)
    stats <- stats[ord]
    statsAdj <- sign(stats) * (abs(stats)^gseaParam)
    statsAdj <- statsAdj / max(abs(statsAdj))
    pathway <- unname(as.vector(na.omit(match(gset, names(statsAdj)))))
    pathway <- sort(pathway)
    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, returnAllExtremes = TRUE)
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops

    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
    diff <- (max(tops) - min(bottoms)) / 8
    es <- c(tops, bottoms)[which.max(abs(c(tops, bottoms)))]

    # run fgsea if no nes and qv given
    if (any(is.na(c(nes, qv)))) {
        gseares <- fgsea(pathways = list(a = gset), stats = stats, nperm = 10000)
        nes <- gseares$NES[1]
        qv <- gseares$padj
    }
    txt <- sprintf("NES : %.2f  \nFDR : %.2e  ", nes, qv)

    # set y axis ticks
    if (is.null(ylims)) ylims <- c(floor(min(bottoms) * 20), ceiling(max(tops) * 20)) / 20

    if ((ylims[2] - ylims[1]) > 0.55) {
        br <- 0.2
        ln1 <- seq(ceiling(ylims[1] * 5) / 5, floor(ylims[2] * 5) / 5, br)
    } else {
        br <- 0.1
        ln1 <- seq(ceiling(ylims[1] * 10) / 10, floor(ylims[2] * 10) / 10, br)
    }
    txt_y_pos <- ifelse(es < 0, max(c(ln1, tops)), min(c(ln1, bottoms)))
    max_es <- ifelse(es > 0, max(tops), min(bottoms))

    xlims <- c(0, ceiling(length(stats) / 500) * 500)
    half_line <- base_size / 2

    # running score
    g1 <- ggplot(toPlot, aes(x = x, y = y)) +
        geom_point(color = line.col, size = 0.1) +
        geom_hline(yintercept = 0, colour = "black", linetype = "dashed", size = lwd * 0.8) +
        geom_hline(yintercept = max_es, colour = "red", linetype = "dashed", size = lwd * 0.8) +
        geom_line(colour = line.col, size = lwd) +
        annotate("text", x = max(toPlot$x), y = txt_y_pos, label = txt, hjust = 1, vjust = .5, size = rel(3.0)) +
        scale_y_continuous(breaks = ln1, limits = ylims) +
        xlim(xlims[1], xlims[2]) +
        labs(y = ytitle, title = mtitle) +
        theme_common_gsea(base_size = base_size) +
        theme(
            plot.title = element_text(hjust = 0.5, vjust = 0.2, face = "bold", margin = unit(c(0, 0, 1.5, 0), "mm")),
            panel.grid.major.y = element_line(colour = "grey85", linetype = "dashed", size = lwd * 0.65),
            axis.title.y = element_text(face = "bold", angle = 90, margin = unit(c(0, 1.5, 0, 0), "mm"), size = rel(0.95)),
            axis.text.y = element_text(size = rel(0.8), colour = "black"),
            plot.margin = margin(half_line, half_line * 2.5, 0, half_line)
        )

    # gene set bar
    g2 <- ggplot(data.frame(x = pathway), aes(x = x, y = -diff / 3, xend = x, yend = diff / 3)) +
        geom_segment(size = ticksSize, colour = "grey35") +
        scale_x_continuous(breaks = c(seq(xlims[1], xlims[2], 5000), xlims[2]), limits = xlims) +
        theme_common_gsea(base_size = base_size) +
        theme(plot.margin = margin(0, half_line * 2.5, half_line, half_line))

    # rank score bar (t stat, logFC etc)
    if (statbar) {
        lead <- toPlot[which(toPlot$y == es), ]
        g1 <- g1 + geom_segment(x = lead$x, xend = lead$x, y = 0, yend = es, color = "red", linetype = 3, size = lwd * 0.25)
        g2 <- g2 + theme(plot.margin = margin(0, half_line * 2.5, 0, half_line))
        g3 <- ggplot(data.frame(x = seq_along(stats), stat = stats), aes(x = x, y = 0, xend = x, yend = stat)) +
            geom_segment(size = ticksSize, colour = "grey35") +
            scale_x_continuous(breaks = c(seq(xlims[1], xlims[2], 5000), xlims[2]), limits = xlims) +
            theme_common_gsea(base_size = base_size) +
            theme(plot.margin = margin(0, half_line * 2.5, half_line, half_line))
    }

    gr1 <- ggplotGrob(g1)
    gr2 <- ggplotGrob(g2)

    if (statbar) {
        gr3 <- ggplotGrob(g3)
        ht.ratio <- c(6, 1, 1)
        gr <- rbind(gr1, gr2, gr3)
    } else {
        ht.ratio <- c(7, 1)
        gr <- rbind(gr1, gr2)
    }

    # gr$widths <- grid::unit.pmax(gr1$widths, gr2$widths, gr3$widths)
    gr$widths <- grid::unit.pmax(gr1$widths, gr2$widths)

    # identify the position of the panels within the gtable
    panid <- gr$layout$t[grep(pattern = "panel", gr$layout$name)]
    gr$heights[panid] <- unit(ht.ratio, "null")

    if (draw) {
        grid.newpage()
        grid.draw(gr)
    }

    return(gr)
}


#' @describeIn plotEnrichment2
#' ggplot theme used in plotEnrichment2
#' @export

theme_common_gsea <- function(base_size = 7) {
    half_line <- base_size / 2
    .theme <- theme(
        text = element_text(face = "plain", size = base_size, colour = "black", family = "Arial"),
        plot.title = element_text(size = rel(1.0)),
        axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        complete = TRUE
    )
    .theme
}
