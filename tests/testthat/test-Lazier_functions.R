test_that("removeNAs works", {
  expect_equal(removeNAs(c(1, NA, 3, NA)), c(1, 3))
  expect_equal(removeNAs(c(NA, NA, NA)), numeric(0))
  expect_equal(removeNAs(c(1, 2, 3)), c(1, 2, 3))
})

test_that("removeNAs with non-numeric input", {
  expect_equal(removeNAs(c("a", NA, "b", NA)), c("a", "b"))
  expect_equal(removeNAs(c(NA, NA, NA)), character(0))
  expect_equal(removeNAs(c("x", "y", "z")), c("x", "y", "z"))
})

test_that("countUni works", {
  expect_equal(countUni(c(1, 2, 2, 3)), 3)
  expect_equal(countUni(c("a", "b", "a", "c")), 3)
  expect_equal(countUni(c(NA, NA, NA)), 1)
})