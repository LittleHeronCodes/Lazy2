#' Transparent theme for ggplot
#'
#' Modified transparent theme from ggpubr
#' @param base_size font base size. default 12
#' @param x_angle rotate x labels
#' @export
#' @examples
#' \dontrun{
#' theme_set(theme_transparent2())
#' }

theme_transparent2 <- function(base_size=12, x_angle=0) {
	half_line <- base_size/2
	xjust <- 0
	if(x_angle > 5) xjust <- 1
	.theme <- theme(
		text = element_text(face = 'bold', size = base_size),# family='Sans'),
		panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		panel.border = element_blank(),
		plot.title = element_text(hjust = 0.5, vjust = 1.5, size = rel(1), margin = unit(c(0.5,0,0.5,0), 'mm')),
		axis.line = element_line(colour = 'black', size = 1),
		axis.text.x = element_text(size = rel(1), colour = 'black'),
		axis.text.y = element_text(size = rel(1), colour = 'black', face = 'plain'),
		legend.key = element_blank(), 
		panel.background = element_rect(fill = "transparent", colour = NA), panel.ontop = TRUE,
		plot.background = element_rect(fill = "transparent",colour = NA),
		legend.background = element_rect(fill = "transparent", colour = NA), 
		strip.background = element_rect(fill = "#F2F2F2", colour = "black", size = 0.7), 
		# plot.margin = margin(half_line, half_line, half_line, half_line), # to reduce space
		plot.margin = unit(c(1,1,1,1), 'mm'),
		complete = TRUE)
	if(x_angle != 0)
		.theme <- .theme + theme(axis.text.x = element_text(angle = x_angle, hjust = xjust, vjust = xjust))
	.theme
}

