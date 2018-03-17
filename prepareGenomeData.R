# prepareGenomeData.R
#
# Purpose:  Read source data from provider databases, map to HGNC symbols
#           and store as tab-delimied text.
#
#           Currently: linear data from HGNC only
#
# Version:  1.0
# Date:     2018 03 04
# Author:   Boris Steipe <boris.steipe@utoronto.ca>
#
# Dependencies:
#           packages: readr, stringr
#
# License: GPL-3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Version history:
#
#   1.0     ADD HPA (Human Protein Atlas) data; add comments on data sources
#   0.5     Add GWAS data
#   0.4     Major effort to properly parse GO and annotate all Chr 20
#              genes with the most informative GOslim term
#   0.3.1   Add gene type to basic data. Use HGNC symbols as authoritative
#           source for gene annotations
#   0.3     Remove BioMart GO data from basic gene data, prepare separate
#             GO data file and populate Gene data with GOslim from BP
#             ontology only.
#   0.2     add STRING data
#   0.1     First draft: Gene data from BioMart
#
# ToDo:
#    - add GWAS and Protein data
#
# ==============================================================================


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                                               Line
#TOC> -------------------------------------------------------------------
#TOC>   1        INIT                                                  62
#TOC>   1.1        Parameters                                          65
#TOC>   1.2        Packages                                            80
#TOC>   1.3        Functions                                           94
#TOC>   2        HGNC SYMBOLS AND CROSSREFERENCES                     123
#TOC>   3        BIOMART GENE ANNOTATIONS                             178
#TOC>   4        STRING DATA                                          224
#TOC>   5        GO DATA                                              291
#TOC>   5.1        GO annotations (for Chr 20 genes)                  312
#TOC>   5.2        Analyze the GO graph                               367
#TOC>   5.2.1          Parse GO term definitions and edges            380
#TOC>   5.2.2          Compile annotation counts                      462
#TOC>   5.3        Fetch GOslim terms                                 558
#TOC>   5.4        Annotate Chr 20 Genes with unique GO terms         576
#TOC>   6        GWAS (GENOME WIDE ASSOCIATION STUDIES)               633
#TOC>   7        HUMAN PROTEIN ATLAS DATA                             676
#TOC>   8        FINISH                                               719
#TOC> 
#TOC> ==========================================================================


# =    1  INIT  ================================================================


# ==   1.1  Parameters  ========================================================

# Paths to the directories that contains the various data sets. The source files
# are described below, but some of them are quite large so we are not
# loading them into the repo. Download information is in the respective
# script sections

HGNCDIR    <- "./"   # Human Gene Nomenclature Committe
BIOMARTDIR <- "./"   # Ensembl Biomart gene annotation
STRINGDIR  <- "./"   # STRING functional interaction data
GODIR      <- "./"   # Gene Ontology data
GWASDIR    <- "./"   # Genome Wide Association Studies
HPADIR     <- "./"   # Human Protein Atlas


# ==   1.2  Packages  ==========================================================
# Load all required packages.

if (!require(readr, quietly=TRUE)) {
  install.packages("readr")
  library(readr)
}

if (!require(stringr, quietly=TRUE)) {
  install.packages("stringr")
  library(stringr)
}


# ==   1.3  Functions  =========================================================

# Utility functions ...

# Some of the data sets we process are quite large. A progress bar tells us
# how we are doing ...
pBar <- function(i, l, nCh = 50) {
  # Draw a progress bar in the console
  # i: the current iteration
  # l: the total number of iterations
  # nCh: width of the progress bar
  ticks <- round(seq(1, l-1, length.out = nCh))
  if (i < l) {
    if (any(i == ticks)) {
      p <- which(i == ticks)[1]  # use only first, in case there are ties
      p1 <- paste(rep("#", p), collapse = "")
      p2 <- paste(rep("-", nCh - p), collapse = "")
      cat(sprintf("\r|%s%s|", p1, p2))
      flush.console()
    }
  }
  else { # done
    cat("\n")
  }
}




# =    2  HGNC SYMBOLS AND CROSSREFERENCES  ====================================

# The HGNC (Human Gene Nomeclature Committee) is the authoritative source for
# recognized genes. Source data is downloaded from the custom download page
# at https://www.genenames.org/cgi-bin/download
#
# - HGNC Approved Symbol
# - Approved Name
# - Locus Type
# - Chromosome
# - Ensembl Gene ID
# - Entrez Gene ID (external)
# - OMIM ID (external)
# - RefSeq ID (external)
# - UniProt ID (external)
# - Ensembl ID (external)
# - UCSC ID (external)

# Download only "Approved" status, and Chr 20

# read source data
tmp <- read_tsv(paste0(HGNCDIR, "HGNC_data.tsv"))   # 1,041 rows

# check which Locus types appear in this data
unique(tmp$`Locus Type`)

# subset rows with gene types of interest
myLocusTypes <- c("gene with protein product",
                  "RNA, transfer")

tmp <- tmp[tmp$`Locus Type` %in% myLocusTypes, ]    # 529 rows

# subset columns of interest
Chr20GeneData <- data.frame(sym = tmp$`Approved Symbol`,
                            name = tmp$`Approved Name`,
                            type = tmp$`Locus Type`,
                            EnsemblID = tmp$`Ensembl Gene ID`,
                            EntrezGeneID = tmp$`Entrez Gene ID(supplied by NCBI)`,
                            OMIMID = tmp$`OMIM ID(supplied by OMIM)`,
                            RefSeqID = tmp$`RefSeq(supplied by NCBI)`,
                            UniProtID = tmp$`UniProt ID(supplied by UniProt)`,
                            UCSCID = tmp$`UCSC ID(supplied by UCSC)`,
                            stringsAsFactors = FALSE)

rownames(Chr20GeneData) <- Chr20GeneData$sym

Chr20GeneData$type <- gsub("gene with protein product",
                           "protein",
                           Chr20GeneData$type)
Chr20GeneData$type <- gsub("RNA, transfer",
                           "tRNA",
                           Chr20GeneData$type)



# =    3  BIOMART GENE ANNOTATIONS  ============================================

# The Ensembl Biomart system provides extensive annotations for genome-sequenced
# model organisms.

# Read file obtained via custom download from ensembl biomart
# http://useast.ensembl.org/
#
# Click on biomart and select:
#
# Ensembl Genes 91
#
# Dataset: Human genes (GRCh38.p10)
#
# Filters: REGION: Chromosome/scaffold: 20
# Attrib.: Gene stable ID
#          Gene start (bp)
#          Gene end (bp)
#          Strand
#          HGNC symbol

# read source data
tmp <- read_tsv(paste0(BIOMARTDIR, "mart_export.txt"))   # 5,379 rows

# use only rows that match Ensembl IDs in Chr20GeneData
tmp <- tmp[(tmp$`Gene stable ID` %in% Chr20GeneData$EnsemblID), ]  # 521 rows

# Add gene start, end, and strand information:
# Initialize columns
Chr20GeneData$start <- as.numeric(NA)
Chr20GeneData$end <- as.numeric(NA)
Chr20GeneData$strand <- as.numeric(NA)

# Add data
for (i in 1:nrow(Chr20GeneData)) {
  iRow <- which(tmp$`Gene stable ID` == Chr20GeneData$EnsemblID[i])
  if (length(iRow) == 1) {
    Chr20GeneData$start[i] <- tmp$`Gene start (bp)`[iRow]
    Chr20GeneData$end[i] <- tmp$`Gene end (bp)`[iRow]
    Chr20GeneData$strand[i] <- tmp$Strand[iRow]
  }
}




# =    4  STRING DATA  =========================================================

# STRING is a database of functional interactions. Interactions are scored
# and made available as network edges.

# Source data was downloaded from STRING database via organism specific
# download.
#
# https://string-db.org/
#
#   9606.protein.links.v10.5.txt  (522.2 MB)  - contains relationships
#   9606.protein.aliases.v10.5.txt (157.7 MB) - contains IDs
#
#

# Read the alias data: we need that to find which ENSP IDs map to which
# HGNC symbols
tmp <- read_tsv(paste0(STRINGDIR, "9606.protein.aliases.v10.5.txt"),
                skip = 1,
                col_names = c("ENSP", "ID", "source"))  # 2,055,779 rows

tmp <- tmp[grep("BioMart_HUGO", tmp$source), ]  # 19,119 rows

tmp <- tmp[tmp$ID %in% Chr20GeneData$sym, ] # 497 of 521 symbols mapped
ENS2symMap <- tmp$ID                        # extract symbols...
names(ENS2symMap) <- tmp$ENSP               # ... and use ENSP IDs as names


# Read the interaction graph data: this is a weighted graph defined as an
# edge list with gene a, gene b, confidence score (0, 999).
tmp <- read_delim(paste0(STRINGDIR, "9606.protein.links.v10.5.txt"),
                  delim = " ",
                  skip = 1,
                  col_names = c("a", "b", "score"))  # 11,353,056 rows

tmp <- tmp[tmp$score >= 900, ]  # 547,620 rows of high-confidence edges

# Extract edges where both genes are in the mapped Chr 20 genes
Chr20FuncIntx <- tmp[(tmp$a %in% names(ENS2symMap)) &
                     (tmp$b %in% names(ENS2symMap)), ]  # 564 rows

# Use ENS2symMap to translate ENSP IDs to HGNC symbols
Chr20FuncIntx$a <- ENS2symMap[Chr20funcIntx$a]
Chr20FuncIntx$b <- ENS2symMap[Chr20funcIntx$b]

# We treat this as an undirected graph, thus we remove duplicates.

# Add a sorted key column

Chr20FuncIntx$keys <- character(nrow(Chr20FuncIntx))
for (i in 1:nrow(Chr20FuncIntx)) {
  Chr20FuncIntx$keys[i] <- paste(sort(c(Chr20FuncIntx$a[i],
                                        Chr20FuncIntx$b[i])),
                                 collapse = ":")
}

# remove rows with duplicated keys
Chr20FuncIntx <- Chr20FuncIntx[! duplicated(Chr20FuncIntx$keys), ]

# remove unneeded column
Chr20FuncIntx <- Chr20FuncIntx[ , c("a", "b", "score")]

# Done
write_tsv(Chr20FuncIntx, path = "Chr20FuncIntx.tsv")



# =    5  GO DATA  =============================================================

# GO (Gene Ontology) is a very large project that has undertaken to organize
# biomolecular concepts into three separate ontologies (Cellular Component,
# Molecular Function, and Biological Process), as well as annotating model
# organism genes to those terms. GO data is large, complex, and has subtle
# pitfalls, and parsing GO terms takes by far the most preprocessing in order to
# map the information into simple per-gene annotations.


# We analyze GO annotations to find the most informative
# GOslim term for each gene. We proceed in two steps:
#  - first we fetch all GOA annotations for Chr 20 genes.
#  - Then we analyze the GOslim graph, split it into the three component
#    ontologies, and determine the number of genes annotated to each term
#    and its children. This allows us to pick the most informative term
#    (least number of annotations) for each gene. We add that term to the
#    Gene data table.



# ==   5.1  GO annotations (for Chr 20 genes)  =================================

# Source data is "goa_human.gaf" (74.5 MB) from
# ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/
#


tmp <- read_tsv(paste0(GODIR, "goa_human.gaf"),
                comment = "!",
                col_names = c("DB",
                              "DB_Object_ID",
                              "Symbol",                 # Gene symbol
                              "Qualifier",
                              "GO_ID",                  # GO ID
                              "DB_Reference",
                              "Evidence_Code",
                              "With_(or)_From",
                              "Aspect",                 # GO Ontology
                              "DB_Object_Name",
                              "DB_Object_Synonym",
                              "DB_Object_Type",
                              "Taxon",
                              "Date",
                              "Assigned_By",
                              "Annotation_Extension",
                              "Gene_Product_Form_ID"
                ))  # 440,931 rows

tmp <- tmp[tmp$Taxon == "taxon:9606", ]  # just to make sure (440,395)

# subset symbol, GO ID and "Aspect" columns (Aspect is C, F, or P for cellular
# Component, molecular Function, or biological Process.)
tmp <- tmp[ , c("Symbol", "GO_ID", "Aspect")]

# Many of these are redundant, due to different sources. We assign every
# row a key, and remove duplicates

tmp$key <- sprintf("%s|%s|%s",
                      tmp$Symbol,
                      tmp$GO_ID,
                      tmp$Aspect)

GOannotations <- tmp[! duplicated(tmp$key), c("Symbol",
                                              "GO_ID",
                                              "Aspect")]  # 268,854 rows


# subset terms annotated to symbols in Chr20GeneData
GOdata$Chr20 <- GOannotations[GOannotations$Symbol %in% Chr20GeneData$sym,
                              c("Symbol", "GO_ID", "Aspect")] # 7,353

# how many Chr 20 genes are annotated?
sum(Chr20GeneData$sym %in% Chr20GOdata$Symbol) # 503 of 529


# ==   5.2  Analyze the GO graph  ==============================================

# The GO graph is a DAG. Actually we have three DAGs: molecular function ("F"),
# biological process ("P"), and cellular component ("C"). For gene annotation
# purposes, we would like to group genes according to functional categories. GO
# publishes a subset of terms of particular importance: GOslim terms. There may
# be more than one GOslim term annotated to one gene, thus we would like to be
# able to choose the most informative one. The most informative (most specific)
# term is the one that has the smallest number of other genes annotated to it
# and its descendants. To compile this information, we have to construct the
# entire DAG, record all annotations, and then propagate annotations up the DAG
# to its roots, recording the number of genes annotated to leaves at each step.

# ===   5.2.1  Parse GO term definitions and edges       

# Source data is "go-basic.obo" from
# http://geneontology.org/page/download-ontology  (33.8 MB)


tmp <- readLines(paste0(GODIR, "go-basic.obo"))
iTerms <- which(tmp == "[Term]")  # 47,136 terms
iEnd <- which(tmp == "[Typedef]")[1] # Start of typedef section
iTerms <- c(iTerms, iEnd)            # append to list of term indices

# Initialize a data frame to hold the term definitions
GOdefs <- data.frame(ID = character(),
                     name = character(),
                     namespace = character(),
                     def = character(),
                     stringsAsFactors = FALSE)

# Initialize a data frame to hold the GO graph edges
GOgraph <- data.frame(ID = character(),
                      parentID = character(),
                      namespace = character(),
                      stringsAsFactors = FALSE)

aspects <- c("C", "F", "P")    # map ontology name to single character
names(aspects) <- c("cellular_component",
                    "molecular_function",
                    "biological_process")


N <- length(iTerms) - 1
for (i in 1:N) {       # for each GO term
                       # (not very efficiently written, runs for a few minutes)
  pBar(i, N)           # update progress bar

  # fetch records for this term
  first <- iTerms[i] + 1
  last <- iTerms[i+1] - 1
  term <- paste(first:last], collapse = "|")

  if (grepl("is_obsolete: true", term)) {
    next
  }

  # parse information
  thisID        <- str_match_all(term, "^id: (GO:\\d+)\\|")[[1]][ ,2]
  thisName      <- str_match_all(term, "\\|name: (.+?)\\|")[[1]][ ,2]
  thisNamespace <- str_match_all(term, "\\|namespace: (.+?)\\|")[[1]][ ,2]
  thisNamespace <- as.character(aspects[thisNamespace])  # make single character
  thisDef       <- str_match_all(term, "\\|def: \"(.+?)\"")[[1]][ ,2]
  m <- str_match_all(term, "\\|is_a: (GO:\\d+) ")[[1]]
  if (length(m) > 0) {
    theseParents  <- m[ ,2]
  } else {
    theseParents <- paste(thisNamespace, "root")
  }

  GOdefs <- rbind(GOdefs,   # add term to data frame
                  data.frame(ID = thisID,
                             name = thisName,
                             namespace = thisNamespace,
                             def = thisDef,
                             stringsAsFactors = FALSE))

  nParents <- length(theseParents)      # add GO graph edges to data frame
  GOgraph <- rbind(GOgraph,
                   data.frame(ID = rep(thisID, nParents),
                              parentID = theseParents,
                              namespace = rep(thisNamespace, nParents),
                              stringsAsFactors = FALSE))
}

rownames(GOdefs) <- GOdefs$ID

# Here are the three root IDs
GOgraph$ID[grep("root", GOgraph$parentID)]
# "GO:0003674" "GO:0005575" "GO:0008150"
GOdefs["GO:0003674",]
GOdefs["GO:0005575",]
GOdefs["GO:0008150",]


# ===   5.2.2  Compile annotation counts                 

# GO graphs are DAGs not trees, thus we can't simply
# propagate counts up to the root, we have to keep track of the actual
# annotated terms to avoid double-counting. We store terms in a list, genes in
# a vector for each term, and whenever we add genes from descendants to a term,
# we unique() the vector.

GOgenes <- list()

# Step 1: Compile cumulative gene lists for Go Terms
N <- nrow(GOannotations)
for (i in 1:N) {                                  # for each annotation
  pBar(i, N)
  thisID <- GOannotations$GO_ID[i]                # fetch term ...
  thisGene <- GOannotations$Symbol[i]             # ... and gene
  if (length(GOgenes[[thisID]]$genes) == 0) {     # initialize if new
    GOgenes[[thisID]]$genes <- thisGene
  } else {                                        # else add gene to vector
    GOgenes[[thisID]]$genes <- c(GOgenes[[thisID]]$genes, thisGene)
  }
}

# Make the gene lists unique
N <- length(GOgenes)
for (i in 1:N) {
  pBar(i, N)
  if (length(GOgenes[[i]]$genes) > 0) {
    GOgenes[[i]]$genes <- unique(GOgenes[[i]]$genes)
  }
}

# To compile counts for all nodes and their descendants, we begin at the leaf
# nodes, propagate terms to their parents and remove them from consideration.
# Then we find all new leaf nodes, and iterate until done.

# Step 2: Propagate gene annotations to root nodes

GOdefs$active <- TRUE  # add a column to flag active nodes

nCycles <- 0  # always add a safetynet when running while() loops  :-)
Nterms <- nrow(GOdefs)
while (sum(GOdefs$active) > 0 && nCycles < 100) {
  nCycles <- nCycles + 1
  cat(sprintf("Cycle: %d, %d active terms.\n",
              nCycles,
              sum(GOdefs$active)))

  for (i in 1:Nterms) {
    pBar(i, Nterms)

    # Identify leafs: a leaf is a node that is active, and not parent to
    # other active node(s).

    if (GOdefs$active[i] == TRUE) {

      thisTerm <- GOdefs$ID[i]
      # find all the node's active children

      children <- GOgraph$ID[thisTerm == GOgraph$parentID]
      sel <- which(GOdefs$ID %in% children)
      activeChildren <- GOdefs$ID[which(GOdefs$active[sel])]

      if (length(activeChildren) == 0) { # it's a leaf
        # propagate annotated genes to parents if any exist
        sel <- which(thisTerm == GOgraph$ID)
        theseParents <- GOgraph$parentID[sel]
        for (parent in theseParents) {
          if (! grepl("root", parent)) {
            GOgenes[[parent]]$genes  <- unique(c(GOgenes[[parent]]$genes,
                                                 GOgenes[[thisTerm]]$genes))
          }
        }
        # unset the "active" flag for this node
        GOdefs[thisTerm, "active"] <- FALSE
      }
    }
  }
}

# Size of GOgenes now: ~ 99 Mb

# Step3: compile the actual counts

GOdefs$counts <- 0
N <- length(GOgenes)
termNames <- names(GOgenes)
for (i in 1:N) {
  pBar(i, N)
  GOdefs[termNames[i], "counts"] <- length(GOgenes[[i]]$genes)
}

# Remove unneeded columns.
GOdefs <- GOdefs[ , c("ID", "name", "namespace", "def", "counts")]


# ==   5.3  Fetch GOslim terms  ================================================

# Source data is "goslim_generic.obo" (257 kb) from
# http://geneontology.org/ontology/subsets/goslim_generic.obo
#

tmp <- readLines(paste0(GODIR, "goslim_generic.obo"))
iTerms <- which(tmp == "[Term]")  # 149 terms

GOslimIDs <- gsub("id: ", "", tmp[iTerms + 1])

# confirm
all(GOslimIDs %in% GOdefs$ID)

# Write GOslim information to file:
write_tsv(GOdefs[GOslimIDs, ], path = "Chr20GOslimData.tsv")


# ==   5.4  Annotate Chr 20 Genes with unique GO terms  ========================


# For each Chr 20 gene, define which GOslim terms it can be annotated to -
# either annotated to the term, or one of its descendants

GOslimChr20 <- data.frame(sym = character(),    # initialize
                          GOslim_ID = character(),
                          namespace = character(),
                          stringsAsFactors = FALSE)

N <- nrow(Chr20GeneData)
for (i in 1:N) {   # add each symbol/term combination
  pBar(i, N)

  sym <- Chr20GeneData$sym[i]

  for(term in GOslimIDs) {
    if (any(GOgenes[[term]]$genes == sym)) {
      GOslimChr20 <- rbind(GOslimChr20,
                           data.frame(sym = sym,
                                      GOslim_ID = term,
                                      namespace = GOdefs[term, "namespace"],
                                      stringsAsFactors = FALSE))
    }
  }
}

# 5,164 rows


# Finally: define most informative GO term of each ontology as the term with the
# smallest number of annotations to it.

Chr20GeneData$GO_C <- NA
Chr20GeneData$GO_F <- NA
Chr20GeneData$GO_P <- NA

for (thisNS in c("C", "F", "P")) {
  thisCol <- paste0("GO_", thisNS)

  for (i in 1:nrow(Chr20GeneData)) {
    sel <- which(GOslimChr20$sym == Chr20GeneData$sym[i]  &
                   GOslimChr20$namespace == thisNS)
    if (length(sel) == 1) {
      Chr20GeneData[i, thisCol] <- GOslimChr20$GOslim_ID[sel]
    } else if (length(sel) > 1) {
      # More than one GO term annotated in this ontology. Choose the
      # one with the smallest number of gene counts (most specific).
      xID <- GOslimChr20$GOslim_ID[sel]
      xCounts <- GOdefs$counts[GOdefs$ID %in% xID]
      bestTerm <- (xID[xCounts == min(xCounts)])[1] # choose the first element
      Chr20GeneData[i, thisCol] <- bestTerm         #     in case of ties
    }
  }
}

# =    6  GWAS (GENOME WIDE ASSOCIATION STUDIES)  ==============================

# The GWAS catalogue provides a compilation of results from studies that have
# discovered significant genotype/phenotype associations in the human genome.

# Source data is downloaded from the GWAS catalogue
#  by selecting for Chr 20 associations https://www.ebi.ac.uk/gwas

tmp <- read_tsv(paste0(GWASDIR,
                       "gwas-2018-03-02-chr20_1-64444167.tsv"))   # 1,520 rows


# retain only rows with correct Chr ID
tmp <- tmp[tmp$CHR_ID == "20", ] # 1,507 rows

x <- unique(tmp$`DISEASE/TRAIT`)

# Filter traits where the reported CHR_POS falls between start and end of a
# Chr 20 gene

Chr20GWAStraits <- character()
N <- nrow(Chr20GeneData)
for (i in 1:N) {
  pBar(i, N)
  iTraits <- which(tmp$CHR_POS >= Chr20GeneData$start[i] &
                  tmp$CHR_POS <= Chr20GeneData$end[i])
  if (length(iTraits > 0)) {
    traits <- tmp$`DISEASE/TRAIT`[iTraits]
  } else {
    traits <- "-"
  }
  Chr20GWAStraits <- c(Chr20GWAStraits, sprintf("%s\t%s",
                                                Chr20GeneData$sym[i],
                                                traits))
}

Chr20GWAStraits <- unique(Chr20GWAStraits)   # remove duplicates

Chr20GWAStraits <- c("sym\ttrait", Chr20GWAStraits)  # add header

writeLines(Chr20GWAStraits, con = "Chr20GWAStraits.tsv")


# =    7  HUMAN PROTEIN ATLAS DATA  ============================================

# The Human Protein Atlas project is a long running program that annotates gene
# products with their subcellular location, prognostic value, importance for
# cancer, and tissue-specificity.

# Source data was downloaded by selecting chromosome 20 as the search term at
# protein atlas.org
#
# https://www.proteinatlas.org/search/chromosome:20 ()
tmp <- read_tsv(paste0(HPADIR, "chromosome_20.tsv"))  # 531 rows

sum(! (tmp$Gene %in% Chr20GeneData$sym)) # 17 not identical to current HGNC

tmp <- as.data.frame(tmp[tmp$Gene %in% Chr20GeneData$sym, ])
rownames(tmp) <- tmp$Gene

# Subset to columns of greatest interest
tmp <- tmp[ , c("Protein class",
                "Subcellular location",
                "Prognostic p-value",
                "RNA cancer category",
                "RNA tissue category",
                "RNA TS TPM",
                "TPM max in non-specific")]

# add rows that are not in Chr20GeneData with all NA values
missingSym <- Chr20GeneData$sym[! Chr20GeneData$sym %in% row.names(tmp)]
tmp[missingSym, ] <- rep(NA, ncol(tmp))

# add data to Chr20GeneData

sym <- rownames(Chr20GeneData)

Chr20GeneData$HPAclass           <- tmp[sym, "Protein class"]
Chr20GeneData$HPAlocation        <- tmp[sym, "Subcellular location"]
Chr20GeneData$HPAprognostic      <- tmp[sym, "Prognostic p-value"]
Chr20GeneData$HPAcancerCat       <- tmp[sym, "RNA cancer category"]
Chr20GeneData$HPAtissueCat       <- tmp[sym, "RNA tissue category"]
Chr20GeneData$HPAspecificExpr    <- tmp[sym, "RNA TS TPM"]
Chr20GeneData$HPAnonSpecificExpr <- tmp[sym, "TPM max in non-specific"]


# =    8  FINISH  ==============================================================

# Final sanity check and cleanup: it turns out some genes do not have
# gene start and end annotated: remove

sel <- which(is.na(Chr20GeneData$start) |
             is.na(Chr20GeneData$end) |
             is.na(Chr20GeneData$strand))

Chr20GeneData <- Chr20GeneData[ -sel, ]

# Gene data annotations completed ... write data frame to file:

write_tsv(Chr20GeneData, path = "Chr20GeneData.tsv")


# [END]
