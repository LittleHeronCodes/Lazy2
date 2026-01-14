
# Plot function for GO enrichment results
draw_plot_ego <- function(pldf, mtitle = "") {
	if (max(pldf$Count) > 10) {
		size_breaks <- seq(5, max(pldf$Count), 5)
	} else {
		size_breaks <- seq(1, max(pldf$Count), 2)
	}

	pldf |> 
		ggplot(aes(x = reorder(Description, FE), y = FE)) +
		geom_point(
			aes(fill = -log10(p.adjust), size = Count), 
			shape = 21, color = "black", stroke = 0.5
		) +
		geom_segment(
			aes(xend = Description, yend = 0), color = "grey", linetype = "dashed"
		) +
		geom_text(
			aes(label = paste0("FE=", round(FE, 2))), hjust = 1.5, size = 3
		) +
		scale_fill_gradientn(
			name = "FDR",
			colours = c("white", "pink", "red", "darkred"),
			breaks = -log10(c(0.2, 0.05, 0.01, 0.001)),
			labels = c("0.2", "0.05", "1e-2", "1e-3"),
			limits = c(0, 3),
			oob = scales::squish
		) +
		scale_size_continuous(
			name = "Count",
			breaks = size_breaks,
			limits = c(1, max(pldf$Count)),
			labels = size_breaks
		) +
		labs(x = "GO term", y = "FE", title = mtitle) +
		coord_flip() +
		theme_minimal() +
		theme(
			legend.position = "bottom",
			axis.text.y = element_text(size = 10)
		) +
		guides(
			fill = guide_colorbar(label.theme = element_text(angle = 90, hjust = 0.5)),
			size = guide_legend()
		)
}
