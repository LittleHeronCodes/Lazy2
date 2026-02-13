##  Get GO gene set information used for clusterProfiler (mild hack)

# libraries
# library(data.table)
library(tidyverse)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Load and parse gson object
##   Function refactored from related functions in clusterProfiler gson_GO

get_GO_terms2 <- function(
	orgdb, keytype = "ENTREZID", ont = "ALL", set_sizes = c(10, 500)
) {
	require(AnnotationDbi)
	require(GO.db)

	if (is(orgdb, "character")) {
		require(orgdb, character.only = TRUE)
		orgdb <- eval(parse(text = orgdb))
	}

	if (is.null(set_sizes) || length(set_sizes) != 2) {
		set_sizes <- c(0, Inf)
	}

	if (!ont %in% c("BP", "MF", "CC", "ALL")) {
		stop("ont should be one of 'BP', 'MF', 'CC', or 'ALL'")
	} 
	
	goterms <- AnnotationDbi::Ontology(GO.db::GOTERM) # all GO terms, with BP/MF/CC
	if (ont != "ALL") {
		goterms <- goterms[goterms == ont]
	}
	
	# get GO terms
	go2gene <- suppressMessages(
		AnnotationDbi::mapIds(
			orgdb, keys = names(goterms), column = keytype, 
			keytype = "GOALL", multiVals = "list"
		)
	)

	# GO term gene list
	gsid2gene <- lapply(go2gene, \(x) unique(x[!is.na(x)]))
	gsize <- sapply(gsid2gene, length)
	gsid2gene <- gsid2gene[which(gsize >= set_sizes[1] & gsize <= set_sizes[2])]

	# get GO term names
	termname <- AnnotationDbi::Term(GO.db::GOTERM)
	goterm <- AnnotationDbi::Ontology(GO.db::GOTERM)

	gsid2name <- data.frame(
		gsid = names(gsid2gene), 
		ont = goterm[names(gsid2gene)],
		Description = termname[names(gsid2gene)],
		row.names = NULL
	)

	out <- list(
		go_terms = gsid2name,
		go_glist = gsid2gene
	)
}

system.time(
	goterm_full2 <- get_GO_terms2(org.Hs.eg.db, keytype = "SYMBOL", ont = "ALL")
)

go_glist <- goterm_full$go_glist

gsize <- sapply(go_glist, length)
go_glist <- go_glist[which(gsize >= 10 & gsize <= 500)]




gsize.keep <- which(gsize >= 10 & gsize <= 500)
go_glist.keep <- go_glist[gsize.keep]




egols_rds <- "network_module/go_enrichment_rna_modules_component.rds"
egols_clean <- readRDS(egols_rds)

featspace <- readLines("variance_filtered_features.txt")
gspace <- featspace[!grepl("__", featspace)]


# check
# all(egols_clean$Description %in% goterm_full$Description)

go_glist <- goterm_full |>
	filter(gsid %in% egols_clean$ID) |>
	(\(x) split(x$geneSymbol, x$gsid))()

go_glist <- lapply(go_glist, \(x) intersect(x, gspace))



library(GO.db)

goid <- c("GO:0070821", "GO:1990266", "GO:0101003")
parents <- purrr::map(goid, ~GO.db::GOCCANCESTOR[[.]]) |> setNames(goid)
