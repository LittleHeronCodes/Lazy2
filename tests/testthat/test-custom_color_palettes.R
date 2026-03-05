test_that("custom_color_palette() returns color palette function", {
  expect_true(is.function(custom_color_palette("OceanFire")))
  expect_true(all(grepl("^#", custom_color_palette("OceanFire")(11))))
})


