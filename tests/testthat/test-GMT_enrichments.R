# Tests for GMT_enrichments.R

test_that("readGMT reads GMT file as list", {
  # Create temporary GMT file
  temp_gmt <- tempfile(fileext = ".gmt")
  writeLines(c(
    "PATHWAY1\tDescription1\tGeneA\tGeneB\tGeneC",
    "PATHWAY2\tDescription2\tGeneD\tGeneE\tGeneF\tGeneG"
  ), temp_gmt)
  
  # Read GMT file
  result <- readGMT(temp_gmt, as.df = FALSE)
  
  # Check structure
  expect_type(result, "list")
  expect_length(result, 2)
  expect_named(result, c("PATHWAY1", "PATHWAY2"))
  
  # Check content
  expect_equal(result$PATHWAY1, c("GeneA", "GeneB", "GeneC"))
  expect_equal(result$PATHWAY2, c("GeneD", "GeneE", "GeneF", "GeneG"))
  
  # Clean up
  unlink(temp_gmt)
})

test_that("readGMT reads GMT file as data frame", {
  # Create temporary GMT file
  temp_gmt <- tempfile(fileext = ".gmt")
  writeLines(c(
    "PATHWAY1\tDescription1\tGeneA\tGeneB",
    "PATHWAY2\tDescription2\tGeneC\tGeneD\tGeneE"
  ), temp_gmt)
  
  # Read GMT file as data frame
  result <- readGMT(temp_gmt, as.df = TRUE)
  
  # Check structure
  expect_s3_class(result, "data.frame")
  expect_named(result, c("ont", "gene"))
  expect_equal(nrow(result), 5)
  
  # Check content (ont is a factor)
  expect_equal(as.character(result$ont), c("PATHWAY1", "PATHWAY1", "PATHWAY2", "PATHWAY2", "PATHWAY2"))
  expect_equal(as.character(result$gene), c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE"))
  
  # Clean up
  unlink(temp_gmt)
})

test_that("writeGMT creates valid GMT file", {
  # Create test genelist
  genelist <- list(
    PATHWAY1 = c("GeneA", "GeneB", "GeneC"),
    PATHWAY2 = c("GeneD", "GeneE")
  )
  
  # Write to temporary file
  temp_gmt <- tempfile(fileext = ".gmt")
  writeGMT(temp_gmt, genelist)
  
  # Read back the file
  lines <- readLines(temp_gmt)
  
  # Check structure
  expect_length(lines, 2)
  expect_match(lines[1], "^PATHWAY1\t\tGeneA\tGeneB\tGeneC$")
  expect_match(lines[2], "^PATHWAY2\t\tGeneD\tGeneE$")
  
  # Clean up
  unlink(temp_gmt)
})

test_that("writeGMT handles geneset descriptions", {
  genelist <- list(
    PATHWAY1 = c("GeneA", "GeneB"),
    PATHWAY2 = c("GeneC", "GeneD")
  )
  
  # Test with named descriptions
  temp_gmt <- tempfile(fileext = ".gmt")
  descriptions <- c(PATHWAY1 = "Desc1", PATHWAY2 = "Desc2")
  writeGMT(temp_gmt, genelist, descriptions)
  
  lines <- readLines(temp_gmt)
  expect_match(lines[1], "^PATHWAY1\tDesc1\tGeneA\tGeneB$")
  expect_match(lines[2], "^PATHWAY2\tDesc2\tGeneC\tGeneD$")
  
  unlink(temp_gmt)
  
  # Test with single description
  temp_gmt2 <- tempfile(fileext = ".gmt")
  writeGMT(temp_gmt2, genelist, "Single description")
  
  lines2 <- readLines(temp_gmt2)
  expect_match(lines2[1], "^PATHWAY1\tSingle description\tGeneA\tGeneB$")
  expect_match(lines2[2], "^PATHWAY2\tSingle description\tGeneC\tGeneD$")
  
  unlink(temp_gmt2)
})

test_that("writeGMT errors on invalid genelist format", {
  # Test with non-list input
  expect_error(
    writeGMT(tempfile(), "not a list"),
    "Check genelist format"
  )
  
  # Test with list containing non-character vectors
  expect_error(
    writeGMT(tempfile(), list(pathway = c(1, 2, 3))),
    "Check genelist format"
  )
  
  # Test with nested list
  expect_error(
    writeGMT(tempfile(), list(pathway = list(c("A", "B")))),
    "Check genelist format"
  )
})

test_that("writeGMT and readGMT round-trip works", {
  # Create test genelist
  genelist <- list(
    PATHWAY1 = c("GeneA", "GeneB", "GeneC"),
    PATHWAY2 = c("GeneD", "GeneE", "GeneF")
  )
  
  # Write and read back
  temp_gmt <- tempfile(fileext = ".gmt")
  writeGMT(temp_gmt, genelist, "Test description")
  result <- readGMT(temp_gmt, as.df = FALSE)
  
  # Check that gene sets match
  expect_equal(result$PATHWAY1, genelist$PATHWAY1)
  expect_equal(result$PATHWAY2, genelist$PATHWAY2)
  
  # Clean up
  unlink(temp_gmt)
})

test_that("hypergeoTestForGeneset performs basic enrichment test", {
  # Create test data
  query <- c("A", "B", "C")
  refGMT <- list(
    SET1 = LETTERS[1:10],
    SET2 = LETTERS[4:25]
  )
  gspace <- LETTERS
  
  # Run test
  result <- hypergeoTestForGeneset(query, refGMT, gspace)
  
  # Check structure
  expect_s3_class(result, "data.frame")
  expect_named(result, c("ID", "pVal", "logP", "qVal", "logQ", "fe", "tan", "int", "gsRatio", "bgRatio", "intGenes"))
  expect_equal(nrow(result), 2)
  
  # Check ID column
  expect_equal(result$ID, c("SET1", "SET2"))
  
  # Check intersection counts
  expect_equal(result$int[1], 3)  # A, B, C are all in SET1
  expect_equal(result$int[2], 0)  # None of A, B, C are in SET2 (which starts at D)
  
  # Check ratio formats
  expect_equal(result$gsRatio[1], "3/3")
  expect_equal(result$bgRatio[1], "10/26")
})

test_that("hypergeoTestForGeneset filters small gene sets", {
  query <- c("A", "B", "C")
  refGMT <- list(
    LARGE_SET = LETTERS[1:15],
    SMALL_SET = c("A", "B"),
    MEDIUM_SET = LETTERS[1:10]
  )
  gspace <- LETTERS
  
  # Should warn and exclude SMALL_SET
  expect_warning(
    result <- hypergeoTestForGeneset(query, refGMT, gspace, minGeneSet = 10),
    "less than.*10 genes and were excluded"
  )
  
  # Check that only LARGE_SET and MEDIUM_SET remain
  expect_equal(nrow(result), 2)
  expect_true("LARGE_SET" %in% result$ID)
  expect_true("MEDIUM_SET" %in% result$ID)
  expect_false("SMALL_SET" %in% result$ID)
})

test_that("hypergeoTestForGeneset warns when query outside gene space", {
  query <- c("A", "B", "C", "extra1", "extra2")
  refGMT <- list(SET1 = LETTERS[1:10])
  gspace <- LETTERS
  
  # Should warn about items outside gene space
  expect_warning(
    result <- hypergeoTestForGeneset(query, refGMT, gspace),
    "query items were found outside of background space"
  )
  
  # Should still produce results with filtered query
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) >= 1)
})

test_that("hypergeoTestForGeneset errors on empty query", {
  query <- character(0)
  refGMT <- list(SET1 = LETTERS[1:10])
  gspace <- LETTERS
  
  # Empty query should error
  expect_error(
    hypergeoTestForGeneset(query, refGMT, gspace),
    "Query length is zero"
  )
})

test_that("hypergeoTestForGeneset errors when all gene sets filtered", {
  query <- c("A", "B", "C")
  refGMT <- list(
    SET1 = c("A", "B"),
    SET2 = c("C", "D", "E")
  )
  gspace <- LETTERS
  
  # All sets have less than 10 genes
  expect_error(
    suppressWarnings(hypergeoTestForGeneset(query, refGMT, gspace, minGeneSet = 10)),
    "Length of refGMT after filtering is zero"
  )
})

test_that("hypergeoTestForGeneset calculates fold enrichment correctly", {
  query <- LETTERS[1:5]
  refGMT <- list(
    ENRICHED = LETTERS[1:10],    # 5/5 overlap, 10/26 in background
    DEPLETED = LETTERS[15:26]    # 0/5 overlap, 12/26 in background
  )
  gspace <- LETTERS
  
  result <- hypergeoTestForGeneset(query, refGMT, gspace, ef.psc = 0)
  
  # Check fold enrichment for enriched set
  expected_fe_enriched <- (5 / ((10/26) * 5))
  expect_equal(result$fe[result$ID == "ENRICHED"], expected_fe_enriched, tolerance = 1e-6)
  
  # Check that depleted set has lower fold enrichment
  expect_true(result$fe[result$ID == "DEPLETED"] < result$fe[result$ID == "ENRICHED"])
})

test_that("hypergeoTestForGeneset calculates Tanimoto coefficient correctly", {
  query <- LETTERS[1:5]  # A, B, C, D, E
  refGMT <- list(
    PERFECT = LETTERS[1:5],           # Tanimoto = 5/5 = 1
    PARTIAL = LETTERS[1:15],          # Tanimoto = 5/15 = 0.333...
    NO_OVERLAP = LETTERS[20:26]       # Tanimoto = 0
  )
  gspace <- LETTERS
  
  result <- hypergeoTestForGeneset(query, refGMT, gspace, minGeneSet = 5)
  
  # Check Tanimoto coefficients
  expect_equal(result$tan[result$ID == "PERFECT"], 1.0, tolerance = 1e-6)
  expect_true(result$tan[result$ID == "PARTIAL"] > 0 & result$tan[result$ID == "PARTIAL"] < 1)
  expect_equal(result$tan[result$ID == "NO_OVERLAP"], 0.0, tolerance = 1e-6)
})

test_that("hypergeoTestForGeneset handles pseudocount in fold enrichment", {
  query <- c("A", "B", "C")
  refGMT <- list(SET1 = LETTERS[10:20])
  gspace <- LETTERS
  
  # Test with different pseudocounts
  result_no_psc <- hypergeoTestForGeneset(query, refGMT, gspace, ef.psc = 0)
  result_with_psc <- hypergeoTestForGeneset(query, refGMT, gspace, ef.psc = 1)
  
  # Pseudocount should affect fold enrichment
  expect_false(identical(result_no_psc$fe, result_with_psc$fe))
  expect_type(result_with_psc$fe, "double")
})

test_that("hypergeoTestForGeneset adjusts p-values correctly", {
  query <- LETTERS[1:5]
  refGMT <- list(
    SET1 = LETTERS[1:10],
    SET2 = LETTERS[2:12],
    SET3 = LETTERS[15:25]
  )
  gspace <- LETTERS
  
  result <- hypergeoTestForGeneset(query, refGMT, gspace)
  
  # qVal should be adjusted p-values
  expect_true(all(result$qVal >= result$pVal))
  
  # logP and logQ should be negative log10 of pVal and qVal
  expect_equal(result$logP, -log10(result$pVal), tolerance = 1e-6)
  expect_equal(result$logQ, -log10(result$qVal), tolerance = 1e-6)
})

test_that("hypergeoTestForGeneset handles zero intersection correctly", {
  query <- LETTERS[1:3]
  refGMT <- list(
    OVERLAP = LETTERS[1:10],
    NO_OVERLAP = LETTERS[20:26]
  )
  gspace <- LETTERS
  
  result <- hypergeoTestForGeneset(query, refGMT, gspace, minGeneSet = 5)
  
  # Check that zero intersection is handled
  no_overlap_row <- result[result$ID == "NO_OVERLAP", ]
  expect_equal(no_overlap_row$int, 0)
  expect_equal(no_overlap_row$qVal, 1)
  expect_equal(no_overlap_row$intGenes, "")
})

test_that("hypergeoTestForGeneset returns intersected genes", {
  query <- c("A", "B", "C", "D")
  refGMT <- list(
    SET1 = c("A", "B", "X", "Y", "Z", "W", "V", "U", "T", "S"),
    SET2 = c("C", "D", "E", "F", "G", "H", "I", "J", "K", "L")
  )
  gspace <- LETTERS
  
  result <- hypergeoTestForGeneset(query, refGMT, gspace)
  
  # Check intersected genes exist and contain expected genes
  set1_genes <- strsplit(result$intGenes[result$ID == "SET1"], "/")[[1]]
  set2_genes <- strsplit(result$intGenes[result$ID == "SET2"], "/")[[1]]
  
  expect_true(all(c("A", "B") %in% set1_genes))
  expect_true(all(c("C", "D") %in% set2_genes))
  expect_equal(length(set1_genes), 2)
  expect_equal(length(set2_genes), 2)
})

test_that("hypergeoTestForGeneset handles genes not in gene space in refGMT", {
  query <- c("A", "B", "C")
  refGMT <- list(
    SET1 = c("A", "B", "extra1", "extra2", "D", "E", "F", "G", "H", "I", "J", "K")
  )
  gspace <- LETTERS
  
  # Should filter refGMT genes to only those in gspace
  result <- hypergeoTestForGeneset(query, refGMT, gspace)
  
  # Check that result is produced
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) >= 1)
  
  # Background ratio should reflect genes in gspace only (10 genes from SET1 in LETTERS)
  expect_match(result$bgRatio[result$ID == "SET1"], "^10/26$")
})
