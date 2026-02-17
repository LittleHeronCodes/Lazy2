#' Custom Color Palettes
#'
#' Custom color palettes for data visualization.
#' @param palette Name of the color palette to use. Use "list" to see all available palettes or "aliases" to see palette aliases.
#' @return A function that generates a color ramp palette.
#' @export
#' @examples
#' colpal <- custom_color_palette("OceanFire")(11)
#' colpal

custom_color_palette <- function(palette = "OceanFire") {
	pal_list <- list(
		# Blue to Red/Orange palettes
		OceanFire  = c( # was BlRd
			"#0045ff", "#4172ff", "#5e9bff", "#7ac2ff", "#a4e7ff", "#ffffff", 
			"#ffd5c3", "#ffa889", "#f57b56", "#e14d2a", "#c90000"
		),
		SapphireFlame = c( # was BlRd1
			"#0033aa", "#1f5ec8", "#3a89e0", "#5bb3f2", "#8edaff", "#ffffff", 
			"#ffddc8", "#ffb18a", "#ff8562", "#ff4f3b", "#d10023"
		),
		CoolBeach = c( # was BlRd2
			'#1a58cb', '#327dd8', '#42a2e3', '#56c7ed', '#88eaef', '#ffffcf', 
			'#fad4a6', '#f2a98a', '#e67b76', '#d54a68', '#bc0059'
		),
		DuskBeach = c( # was BlRd3
			'#373da8', '#3965c0', '#3f8bd5', '#5ab1e4', '#94d4e8', '#f6edd0', 
			'#f6c4ac', '#ed9b90', '#de7378', '#cc4863', '#b80049'
		),
		CoolAmber = c( # was BlOr
			'#007cff', '#3b98f9', '#69b0f5', '#94c6f2', '#c0daef', '#ededed', 
			'#ebd5b9', '#e6bc88', '#e0a25c', '#db8732', '#d76700'
		),

		# Teal/Green to Orange/Red palettes
		CoralLagoon  = c( # was GrOr
			"#006063", "#097d7f", "#239a9c", "#45b7b8", "#73d4d4", "#e2e2e2", 
			"#f2c4b4", "#f7a78a", "#f48b66", "#eb7049", "#db5735"
		),
		TropicSunset  = c( # was GrRd
			'#004e6d', '#107288', '#2c98a0', '#50bfb0', '#87e5b4', '#fefda7', 
			'#facc84', '#ee9c65', '#da6c4e', '#bf3d3f', '#9a003a'
		),

		# Green to Pink/Purple palettes
		MossBloom = c( # was GrPn
			"#006d32", "#238f57", "#3fb776", "#5ed19a", "#8ef0c1", "#f0f0f0",
			"#f2d1e2", "#eaa4c9", "#e275ad", "#d44a95", "#c2187b"
		),
		EnchantedGrove = c( # was GrPn1
			'#009e00', '#39ac57', '#5bba88', '#80c6ad', '#abd1c9', '#dbdbdb', 
			'#d7c0e3', '#d3a5e9', '#cf88ec', '#cd69ea', '#d739d3'
		),

		# Blue to Brown palettes (Earth tones)
		DeepEarth = c( # was BlBr
			'#0065af', '#2482bb', '#49a0c2', '#73bcc3', '#a7d6c0', '#ece9ce', 
			'#e1c992', '#ccaa6f', '#b98a53', '#a86a3a', '#974924'
		),
		StormyCanyons = c( # was BlBr2
			'#002e60', '#2e5181', '#5777a0', '#849dbd', '#b6c5d7', '#ededed', 
			'#d0bfae', '#af9378', '#8e684a', '#6c4024', '#4a1900'
		),
		DesertOasis = c( # oasis waters to desert sands (copilot suggestion, modified)
			'#001c73', '#31458e', '#566fa9', '#7f99c3', '#acc5dc', '#eeeeee',
			'#dbc1ab', '#c29579', '#a96a4d', '#8e3e26', '#720000'
		),

		# Blue to Yellow palettes
		TwilightMeadow = c( # was BlYl
			'#4e45da', '#5670e4', '#6b95ec', '#8db8f2', '#bcd8f5', '#f5f5f5', 
			'#ece8b8', '#ddda8c', '#cccd69', '#babf4c', '#a9b134'
		),
		SkyHarvest = c( 
			'#0033cc', '#2e5fd1', '#5a8ae0', '#89b4f0', '#bfe0ff', '#fef9e6', 
			'#f7e6b3', '#edcf80', '#e2b952', '#d5a22b', '#c78d0c'
		),
		CometTail = c(
			'#001a99', '#264db3', '#4d7fcc', '#7aa6e0', '#b3d1f2', '#f2f2e6', 
			'#e6d9b3', '#d9c080', '#ccab52', '#bf942b', '#b3800c'
		),
		CometTailBright = c(
			'#006dff', '#498dff', '#7ba9ff', '#a8c5ff', '#d2e0fe', '#fcfcfc', 
			'#fdedbc', '#fadc88', '#f8ca5d', '#f6b735', '#f4a200'
		),

		# Sequential
		GalaxySpiral = c(
			'#000061', '#11028e', '#320ab6', '#5519d5', '#7b2bec', '#a241f7', 
			'#c857fa', '#e572fd', '#f694ff', '#ffb7ff', '#ffdcff'
		),
		ArcticNight = c(
			'#0a0a23', '#1a1a3a', '#2a2a52', '#3a4a6a', '#4a6a82', '#5a8a9a', 
			'#6aaab2', '#7acaca', '#8aeae2', '#9afffa', '#aaffff'
		)
	)

	# Palette aliases for backwards compatibility and set easy defaults for palette family
	palette_aliases <- list(
		"BlRd"  = "OceanFire",
		"BlRd1" = "SapphireFlame",
		"BlRd2" = "CoolBeach",
		"BlRd3" = "DuskBeach",
		"BlOr"  = "CoolAmber",
		"GrOr"  = "CoralLagoon",
		"GrRd"  = "TropicSunset",
		"GrPn"  = "MossBloom",
		"GrPn1" = "EnchantedGrove",
		"BlBr"  = "DeepEarth",
		"BlBr2" = "StormyCanyons",
		"BlYl"  = "TwilightMeadow",
		"BlGd"  = "SkyHarvest"
	)

	if (palette %in% names(palette_aliases)) {
		palette <- palette_aliases[[palette]]
	}

	# Return list of available palettes if requested
	if (!(palette %in% names(pal_list))) {
		if (palette == "list") {
			return(names(pal_list))
		} else if (palette == "aliases") {
			return(names(palette_aliases))
		}
		stop("Available palettes: ", paste(names(pal_list), collapse = ", "))
	}

	grDevices::colorRampPalette(pal_list[[palette]])
}

#' @describeIn custom_color_palette Print all available palettes
#' @param n Number of colors to display for each palette
#' @export
#' @examples 
#' show_custom_palettes(11)

show_custom_palettes <- function(n = 11) {
	pals <- custom_color_palette("list")
	par.default <- par(no.readonly = TRUE)
	par(mfrow = c(ceiling(length(pals) / 2), 2), mar = c(2, 2, 2, 2))
	for (pal_name in pals) {
		pal_colors <- custom_color_palette(palette = pal_name)(n)
		barplot(
			rep(1, length(pal_colors)),
			col = pal_colors,
			main = pal_name,
			border = NA,
			axes = FALSE,
			space = 0
		)
	}
	par(par.default)
}
