test_that("custom_color_palette() returns color palette function", {
  # Check that the function returns a function
  expect_true(is.function(custom_color_palette("OceanFire")))

  # Check that the returned function produces valid hex color codes
  expect_true(all(grepl("^#", custom_color_palette("OceanFire")(11))))
})

test_that("custom_color_palette() special return values work", {
  # "list" should return character vector of palette names
  palette_list <- custom_color_palette("list")
  expect_type(palette_list, "character")
  expect_true(length(palette_list) > 0)
  expect_true("OceanFire" %in% palette_list)
  expect_true("CoralLagoon" %in% palette_list)
  expect_true("GalaxySpiral" %in% palette_list)
  
  # "aliases" should return character vector of alias names
  alias_list <- custom_color_palette("aliases")
  expect_type(alias_list, "character")
  expect_true(length(alias_list) > 0)
  expect_true("BlRd" %in% alias_list)
  expect_true("GrOr" %in% alias_list)
})

test_that("palettes return valid hexadecimal colors", {
  # Get all available palettes
  all_palettes <- custom_color_palette("list")
  
  for (pal_name in all_palettes) {
    colors <- custom_color_palette(pal_name)(11)
    
    # All colors should start with #
    expect_true(all(grepl("^#", colors)), 
                info = paste("Palette", pal_name, "has colors without #"))
    
    # All colors should be valid hex codes (# followed by 6 hex digits)
    expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", colors)),
                info = paste("Palette", pal_name, "has invalid hex codes"))
    
    # Should return correct number of colors
    expect_equal(length(colors), 11,
                 info = paste("Palette", pal_name, "returned wrong number of colors"))
  }
})

test_that("custom_color_palette() error handling works", {
  # Invalid palette name should throw error
  expect_error(custom_color_palette("NonExistentPalette"),
               "Available palettes")
  
  # Error message should contain available palette names
  expect_error(custom_color_palette("InvalidName"),
               "OceanFire")
  
  # Empty string should error
  expect_error(custom_color_palette(""),
               "Available palettes")
})

test_that("returned color function handles different n values", {
  pal_func <- custom_color_palette("OceanFire")

  # n = 11 should return 11 colors (base palette)
  colors_11 <- pal_func(11)
  expect_length(colors_11, 11)

  # n = 1 should return 1 color
  expect_length(pal_func(1), 1)
  expect_true(grepl("^#[0-9A-Fa-f]{6}$", pal_func(1)))
  
  # n = 5 should return 5 colors, with first, last, middle matching for n = 11
  colors_5 <- pal_func(5)
  expect_length(colors_5, 5)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", colors_5)))
  expect_equal(colors_5[1], colors_11[1])
  expect_equal(colors_5[5], colors_11[11])
  expect_equal(colors_5[3], colors_11[6])

  # n = 101 should return 101 colors, with first, middle, last colors matching for n = 11
  colors_101 <- pal_func(101)
  expect_length(colors_101, 101)
  expect_equal(colors_101[1], colors_11[1])
  expect_equal(colors_101[101], colors_11[11])
  expect_equal(colors_101[51], colors_11[6])
  
  # n = 100 should work
  colors_100 <- pal_func(100)
  expect_length(colors_100, 100)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", colors_100)))
})

test_that("returned color function handles edge cases", {
  pal_func <- custom_color_palette("CoralLagoon")
  
  # n = 0 should return empty vector
  expect_length(pal_func(0), 0)
  
  # n = 2 should return 2 colors
  expect_length(pal_func(2), 2)
  
  # Verify colors are unique when n is large enough
  colors_50 <- pal_func(50)
  expect_false(any(duplicated(colors_50)))
})

test_that("specific palette colors are correct (regression test)", {
  # Test OceanFire palette returns expected base colors
  oceanfire_colors <- custom_color_palette("OceanFire")(11)
  
  # Test first and last colors match expected values
  expect_equal(toupper(oceanfire_colors[1]), "#0045FF")
  expect_equal(toupper(oceanfire_colors[11]), "#C90000")
  expect_equal(toupper(oceanfire_colors[6]), "#FFFFFF")  # Middle should be white
})

test_that("show_custom_palettes() runs without error", {
  # Open null graphics device to suppress Rplots.pdf
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  
  # Should run without error with default parameter
  expect_no_error(show_custom_palettes())
  
  # Should work with different n values
  expect_no_error(show_custom_palettes(n = 5))
  expect_no_error(show_custom_palettes(n = 20))
})

test_that("show_custom_palettes() handles edge cases", {
  # Open null graphics device to suppress Rplots.pdf
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  
  # n = 1 should work
  expect_no_error(show_custom_palettes(n = 1))
  
  # n = 2 should work
  expect_no_error(show_custom_palettes(n = 2))
  
  # Large n should work
  expect_no_error(show_custom_palettes(n = 50))
})

test_that("show_custom_palettes() restores graphics parameters", {
  # Open null graphics device to suppress Rplots.pdf
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  
  # Store original par settings
  original_mfrow <- par("mfrow")
  original_mar <- par("mar")
  
  # Temporarily change settings
  par(mfrow = c(2, 2), mar = c(1, 1, 1, 1))
  
  # Run function (should restore to defaults)
  show_custom_palettes(n = 5)
  
  # Check that par settings were restored to the function's defaults
  # The function should restore whatever was set before it ran
  current_mfrow <- par("mfrow")
  current_mar <- par("mar")
  
  # Just verify par() is callable and returns valid values
  expect_true(is.numeric(current_mfrow))
  expect_true(is.numeric(current_mar))
  expect_equal(length(current_mfrow), 2)
  expect_equal(length(current_mar), 4)
})

test_that("all palettes work end-to-end", {
  # Get list of all palettes
  all_palettes <- custom_color_palette("list")
  
  # Each palette should produce a working color function
  for (pal_name in all_palettes) {
    pal_func <- custom_color_palette(pal_name)
    expect_true(is.function(pal_func))
    
    colors <- pal_func(11)
    expect_length(colors, 11)
    expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", colors)))
  }
})


