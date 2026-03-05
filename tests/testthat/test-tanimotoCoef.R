test_that("tanimoto_coef() calculates the Tanimoto coefficient correctly", {
  expect_equal(tanimoto_coef(1:5, 4:8), 0.25)
  expect_equal(tanimoto_coef(1:5, 1:5), 1)
  expect_equal(tanimoto_coef(1:5, 6:10), 0)
})

test_that("fold_enrichment() calculates the fold enrichment correctly", {
  expect_equal(fold_enrichment(1:6, 4:8, 1:10), 1)
  expect_equal(fold_enrichment(1:6, 4:8, 1:10, psc = 1), 1)
  expect_equal(fold_enrichment(1:6, 7:8, 1:10), 0)
})
