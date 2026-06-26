test_that("tanimoto_coef() calculates the Tanimoto coefficient correctly", {
  # Basic Tanimoto coefficient calculations
  expect_equal(tanimoto_coef(1:5, 4:8), 0.25)
  expect_equal(tanimoto_coef(1:5, 1:5), 1)
  expect_equal(tanimoto_coef(1:5, 6:10), 0)
})

test_that("tanimoto_coef() works with background space parameter T", {
  # T parameter restricts A and B to only elements in T
  expect_equal(tanimoto_coef(1:10, 5:15, T = 1:10), tanimoto_coef(1:10, 5:10))
  expect_equal(tanimoto_coef(1:5, 3:7, T = 1:10), 3/7)
  
  # T excludes elements from both sets
  expect_equal(tanimoto_coef(1:5, 6:10, T = 3:8), 0)
  expect_equal(tanimoto_coef(1:5, 4:8, T = 4:5), 1)
})

test_that("tanimoto_coef() handles edge cases", {
  # Empty sets
  expect_equal(tanimoto_coef(c(), c()), NaN)
  expect_equal(tanimoto_coef(1:5, c()), 0)
  expect_equal(tanimoto_coef(c(), 1:5), 0)
  
  # Single element sets
  expect_equal(tanimoto_coef(1, 1), 1)
  expect_equal(tanimoto_coef(1, 2), 0)
  
  # Duplicate elements in input (union/intersect handle these)
  expect_equal(tanimoto_coef(c(1, 1, 2, 2), c(2, 2, 3, 3)), 1/3)
})

test_that("fold_enrichment() calculates the fold enrichment correctly", {
  expect_equal(fold_enrichment(1:6, 4:8, 1:10), 1)
  expect_equal(fold_enrichment(1:6, 4:8, 1:10, psc = 1), 1)
  expect_equal(fold_enrichment(1:6, 7:8, 1:10), 0)
})

test_that("fold_enrichment() handles edge cases", {
  # Empty sets should return NA
  expect_true(is.na(fold_enrichment(c(), 1:5, 1:10)))
  expect_true(is.na(fold_enrichment(1:5, c(), 1:10)))
  
  # No overlap between sets
  expect_equal(fold_enrichment(1:3, 7:9, 1:10), 0)
  
  # Complete overlap
  expect_equal(fold_enrichment(1:5, 1:5, 1:10), 2)
  
  # Pseudocount affects low expected counts
  expect_equal(fold_enrichment(1:2, 1:2, 1:100, psc = 1), (2 + 1) / (2/100 * 2 + 1))
  
  # Sets outside of T are intersected with T first
  expect_equal(fold_enrichment(1:15, 5:20, 1:10), fold_enrichment(1:10, 5:10, 1:10))
})

test_that("hypergeo_test() performs hypergeometric test correctly", {
  # Known test case
  result <- hypergeo_test(1:10, 5:15, 1:20)
  
  # Check output structure
  expect_s3_class(result, "data.frame")
  expect_named(result, c("pVal", "oddsRatio", "int", "bg"))
  expect_equal(nrow(result), 1)
  
  # Check values are reasonable
  expect_true(result$pVal >= 0 && result$pVal <= 1)
  expect_true(result$oddsRatio > 0)
  expect_equal(result$int, length(intersect(1:10, 5:15)))
  expect_equal(result$bg, 20)
  
  # Perfect overlap should give significant result
  result_perfect <- hypergeo_test(1:5, 1:5, 1:100)
  expect_lt(result_perfect$pVal, 0.05)
  expect_gt(result_perfect$oddsRatio, 1)
})

test_that("hypergeo_test() handles edge cases", {
  # No overlap
  result <- hypergeo_test(1:5, 6:10, 1:20)
  expect_equal(result$int, 0)
  expect_equal(result$pVal, 1)
  
  # Empty query
  result <- hypergeo_test(c(), 1:10, 1:20)
  expect_equal(result$int, 0)
  
  # Query and reference outside bgspace are intersected first
  result <- hypergeo_test(1:30, 15:40, 1:20)
  expect_equal(result$bg, 20)
  expect_equal(result$int, length(intersect(1:20, 15:20)))
})

test_that("harmean() calculates harmonic mean correctly", {
  # Basic calculation
  expect_equal(harmean(c(1, 2)), 1 / mean(c(1, 0.5)))
  expect_equal(harmean(c(1, 1, 1)), 1)
  expect_equal(harmean(c(2, 4)), 8/3)
  
  # Single value
  expect_equal(harmean(5), 5)
  
  # Typical p-value use case
  pvals <- c(0.01, 0.02, 0.03)
  expect_type(harmean(pvals), "double")
  expect_true(harmean(pvals) > 0 && harmean(pvals) < 1)
})

test_that("harmean() handles NA values", {
  # With NA present
  expect_true(is.na(harmean(c(1, 2, NA))))
  expect_equal(harmean(c(1, 2, NA), na.rm = TRUE), harmean(c(1, 2)))
})

test_that("harmean() handles edge cases", {
  # Zero should cause 0 result (1/mean(Inf, 1) = 1/Inf = 0)
  expect_equal(harmean(c(0, 1)), 0)
  
  # Empty vector
  expect_true(is.nan(harmean(c())))
})

test_that("deprecated functions still work and warn", {
  # tanimotoCoef is deprecated
  expect_warning(result <- tanimotoCoef(1:5, 3:7), "deprecated")
  expect_equal(result, tanimoto_coef(1:5, 3:7))
  
  # getEnrichmentFactor is deprecated
  expect_warning(result <- getEnrichmentFactor(1:6, 4:8, 1:10), "deprecated")
  expect_equal(result, fold_enrichment(1:6, 4:8, 1:10))
  
  # hypergeoTest is deprecated
  expect_warning(result <- hypergeoTest(1:10, 5:15, 1:20), "deprecated")
  expect_equal(result, hypergeo_test(1:10, 5:15, 1:20))
})
