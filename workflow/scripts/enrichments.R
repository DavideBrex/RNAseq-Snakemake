log <- file(snakemake@log[[1]], open= "wt")
sink(log)
sink(log, type = "message")

library(tidyverse)
library(clusterProfiler)
library(ReactomePA)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(openxlsx)

source("workflow/scripts/custom_functions.R")

#------------------------------------------------------------------------------------------
# Read files and parameters
#------------------------------------------------------------------------------------------
DEA.annot <- read.delim(snakemake@input[[1]])

genome       <- snakemake@params[["genome"]]
pvalue       <- as.numeric(snakemake@params[["pvalue"]])
qvalue       <- as.numeric(snakemake@params[["qvalue"]])
fpkm         <- as.numeric(snakemake@wildcards[["fpkm"]])
set_universe <- as.logical(snakemake@params[["set_universe"]])


# Set mouse or human references for the databases
if (genome == "mouse") { 
  kegg.genome <- "mmu"
  db          <- org.Mm.eg.db
} else if (genome == "human") {
  kegg.genome <- "hsa"
  db          <- org.Hs.eg.db
}

#------------------------------------------------------------------------------------------
# Prepare data
#------------------------------------------------------------------------------------------
# Get UP and DOWN-regulated 
UP <- DEA.annot %>% 
  dplyr::filter(DEG == "Upregulated") %>% 
  dplyr::select(Geneid) %>% 
  pull 

DWN <- DEA.annot %>% 
  dplyr::filter(DEG == "Downregulated") %>% 
  dplyr::select(Geneid) %>% 
  pull 

# Get universe of genes. Genes that have been considered for differential expression.
if (set_universe == TRUE) {
  universe <- DEA.annot %>%
    mutate(max_fpkm = purrr::reduce(dplyr::select(., contains("_FPKM")), pmax)) %>%
    dplyr::filter(max_fpkm > 0) %>% 
    dplyr::select(Geneid) %>% 
    pull %>% unique()
} else {
  universe <- DEA.annot %>%
    dplyr::select(Geneid) %>% 
    pull %>% unique()
}
print(paste("Number of genes in universe: ", length(universe)))
#------------------------------------------------------------------------------------------
# Perform enrichments
#------------------------------------------------------------------------------------------
UP.go.bp      <- goEnrichment(UP, genome = genome, ont = "BP", db = db, pvalue = pvalue, qvalue = qvalue, universe = universe)
DWN.go.bp     <- goEnrichment(DWN, genome = genome, ont = "BP", db = db, pvalue = pvalue, qvalue = qvalue, universe = universe)
UP.go.cc      <- goEnrichment(UP, genome = genome, ont = "CC", db = db, pvalue = pvalue, qvalue = qvalue, universe = universe)
DWN.go.cc     <- goEnrichment(DWN, genome = genome, ont = "CC", db = db, pvalue = pvalue, qvalue = qvalue, universe = universe)
UP.go.mf      <- goEnrichment(UP, genome = genome, ont = "MF", db = db, pvalue = pvalue, qvalue = qvalue, universe = universe)
DWN.go.mf     <- goEnrichment(DWN, genome = genome, ont = "MF", db = db, pvalue = pvalue, qvalue = qvalue, universe = universe)

UP.kegg    <- KEGGenrichment(UP, org = kegg.genome, db = db, pvalue = pvalue, qvalue = qvalue, universe = universe)
DWN.kegg   <- KEGGenrichment(DWN, org = kegg.genome, db = db, pvalue = pvalue, qvalue = qvalue, universe = universe)

UP.reactome      <- ReactomeEnrichment(UP, genome = genome, pvalue = pvalue, qvalue = qvalue, universe = universe)
DWN.reactome     <- ReactomeEnrichment(DWN, genome = genome, pvalue = pvalue, qvalue = qvalue, universe = universe)

UP_generalEnrich <- GeneralEnrichment(UP, genome = genome, pvalue = pvalue, qvalue = qvalue, universe = universe)
DWN_generalEnrich <- GeneralEnrichment(DWN, genome = genome, pvalue = pvalue, qvalue = qvalue, universe = universe)

GSEA.hall  <- GSEA_enrichment(DEA.annot, genome = genome, pathways.gmt="Hallmark")
GSEA.c2all <- GSEA_enrichment(DEA.annot, genome = genome , pathways.gmt="C2/M2")
GSEA.c3tft <- GSEA_enrichment(DEA.annot, genome = genome, pathways.gmt = "C3/M3")


#------------------------------------------------------------------------------------------
# save outputs to xls
#------------------------------------------------------------------------------------------

list_of_datasets <- list("GO-BP Up"                 = UP.go.bp@result %>%  dplyr::select(-c("Description"))  %>% as_tibble(rownames = NULL) , 
                         "GO-BP Down"               = DWN.go.bp@result %>%  dplyr::select(-c("Description"))  %>% as_tibble(rownames = NULL) ,
                         "GO-CC Up"                 = UP.go.cc@result %>%  dplyr::select(-c("Description"))  %>% as_tibble(rownames = NULL) ,
                         "GO-CC Down"               = DWN.go.cc@result %>%  dplyr::select(-c("Description"))  %>% as_tibble(rownames = NULL) ,
                         "GO-MF Up"                 = UP.go.mf@result %>%  dplyr::select(-c("Description"))  %>% as_tibble(rownames = NULL) ,
                         "GO-MF Down"               = DWN.go.mf@result %>%  dplyr::select(-c("Description"))  %>% as_tibble(rownames = NULL) ,
                         "KEGG Up"                  = UP.kegg@result %>%  dplyr::select(-c("Description"))  %>% as_tibble(rownames = NULL) ,
                         "KEGG Down"                = DWN.kegg@result %>%  dplyr::select(-c("Description"))  %>% as_tibble(rownames = NULL) ,
                         "Reactome Up"              = UP.reactome@result %>%  dplyr::select(-c("Description"))  %>% as_tibble(rownames = NULL) ,
                         "Reactome Down"            = DWN.reactome@result %>%  dplyr::select(-c("Description"))  %>% as_tibble(rownames = NULL) ,
                         "TF targets (GTRD) Up"     = UP_generalEnrich$`TF targets - GTRD`@result %>%  dplyr::select(-c("Description"))  %>% as_tibble(rownames = NULL) ,
                         "TF targets (GTRD) Down"   = DWN_generalEnrich$`TF targets - GTRD`@result %>%  dplyr::select(-c("Description"))  %>% as_tibble(rownames = NULL) ,
                         "Cell type Up"             = UP_generalEnrich$`Cell Types`@result %>%  dplyr::select(-c("Description"))  %>% as_tibble(rownames = NULL) ,
                         "Cell type Down"           = DWN_generalEnrich$`Cell Types`@result %>%  dplyr::select(-c("Description"))  %>% as_tibble(rownames = NULL) ,
                         "Hallmarks Up - ORA"      =  UP_generalEnrich$Hallmark@result %>%  dplyr::select(-c("Description"))  %>% as_tibble(rownames = NULL) ,
                         "Hallmarks Down - ORA"    =  DWN_generalEnrich$Hallmark@result %>%  dplyr::select(-c("Description"))  %>% as_tibble(rownames = NULL) ,
                         "Hallmarks - GSEA"         = GSEA.hall %>% dplyr::arrange((padj)) %>% dplyr::select(-c("nMoreExtreme", "ES")),
                         "C2all - GSEA"             = GSEA.c2all  %>% dplyr::arrange((padj)) %>% dplyr::select(-c("nMoreExtreme", "ES"))  ,
                         "C3tft - GSEA"             = GSEA.c3tft %>%  dplyr::arrange((padj)) %>% dplyr::select(-c("nMoreExtreme", "ES"))
                         )



write.xlsx(list_of_datasets, file = snakemake@output[["enrichments"]])
