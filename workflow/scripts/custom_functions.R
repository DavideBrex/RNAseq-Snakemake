#----------------------------------------------------------------------------------------------
# Annotation of differential expression table
#----------------------------------------------------------------------------------------------
Annot_DE <- function(df, log2FC = 2, padjust = 0.05, fpkm = 0) {
  require(dplyr)
  require(purrr)
    
    df <- mutate(df, max_fpkm = reduce(select(df, contains("_FPKM")), pmax)) # Add column with max fpkm value between conditions
    
    df$DEG <- "NS"
    df$DEG[which(df$log2FoldChange >= log2FC & df$padj <= padjust & df$max_fpkm >= fpkm)] <- "Upregulated" #Annotate UPregulated genes
    df$DEG[which(df$log2FoldChange <= -log2FC & df$padj <= padjust & df$max_fpkm >= fpkm)] <- "Downregulated" #Annotate DOWNregulated genes
    
    df <- select(df, -max_fpkm) # Remove temporary max fpkm column
  
  return(df)
}


#----------------------------------------------------------------------------------------------
# Functions to transform counts into fpkm and tpm
#----------------------------------------------------------------------------------------------
do_fpkm = function (counts, effective_lengths) {
  exp(log(counts) - log(effective_lengths) - log(sum(counts)) + log(1E9))
}

do_tpm = function (counts, effective_lengths) {
  rate = log(counts) - log(effective_lengths)
  exp(rate - log(sum(exp(rate))) + log(1E6))
}


#----------------------------------------------------------------------------------------------
# Function to downsample count matrix
#----------------------------------------------------------------------------------------------
Down_Sample_Matrix <-
function (expr_mat) 
{
    min_lib_size <- min(colSums(expr_mat))
    down_sample <- function(x) {
        prob <- min_lib_size/sum(x)
        return(unlist(lapply(x, function(y) {
            rbinom(1, y, prob)
        })))
    }
    down_sampled_mat <- apply(expr_mat, 2, down_sample)
    return(down_sampled_mat)
}


#----------------------------------------------------------------------------------------------
# Enrichment analyses functions
#----------------------------------------------------------------------------------------------
#df --> vector of gene symbols of interest
#genome --> "human" or "mouse"
#ont --> "MF", "BP" or "CC"
goEnrichment <- function(
    df, 
    genome ,
    ont      = "BP", 
    db       = org.Mm.eg.db, 
    pvalue   = 0.05, 
    qvalue   = 0.05,
    universe = NULL
) {
  require(clusterProfiler)
  
  if(genome == "human"){
    pathToFile <- paste0('/shares/CIBIO-Storage/scratch_CIBIO/sharedLC/DavideBressan/RNA-seq/ALL_gene_sets/',genome,'/c5.go.',tolower(ont),'.v2024.1.Hs.symbols.gmt')
  } else if(genome == "mouse"){
    pathToFile <- paste0('/shares/CIBIO-Storage/scratch_CIBIO/sharedLC/DavideBressan/RNA-seq/ALL_gene_sets/',genome,'/m5.go.',tolower(ont),'.v2024.1.Mm.symbols.gmt')
  }
  
  stopifnot(file.exists(pathToFile))
  gs_go <- read.gmt(pathToFile)
  colnames(gs_go) <- c("term","gene")
  #now we use trick to fix clusterProfiler background problem https://github.com/YuLab-SMU/clusterProfiler/issues/283
  #trick from https://academic.oup.com/bioinformaticsadvances/article/4/1/vbae159/7829164#493480839
  # we add all genes to the universe as "background" term to avoid the problem
  if(!is.null(universe)){
    bgdf <- data.frame("background",universe)
    colnames(bgdf) <- c("term","gene")
    gs_go <- rbind(gs_go,bgdf)
  }  
  
  ora_res <- enricher(gene = df ,
                      universe = universe, 
                      minGSSize = 5, 
                      maxGSSize = 500000, 
                      TERM2GENE = gs_go,
                      pAdjustMethod="fdr", 
                      pvalueCutoff = 1,
                      qvalueCutoff = 1  )
  
  #add enrichment score
  
  gr <- as.numeric(sapply(strsplit(ora_res@result$GeneRatio,"/"),"[[",1)) /
    as.numeric(sapply(strsplit(ora_res@result$GeneRatio,"/"),"[[",2))
  
  br <- as.numeric(sapply(strsplit(ora_res@result$BgRatio,"/"),"[[",1)) /
    as.numeric(sapply(strsplit(ora_res@result$BgRatio,"/"),"[[",2))
  
  ora_res@result$es <- gr/br
  #ora_res@result$Description=NULL
  
  #now we fix the wrong FDR calculation
  nsets <- length(which(table(gs_go$term)>=5))
  nres <- nrow(ora_res@result)
  if (nres == 0) {
    return(data.frame())
  }
  diff <- nsets - nres
  pvals <- c(ora_res@result$pvalue,rep(1,diff))
  ora_res@result$p.adjust <- p.adjust(pvals,method="fdr")[1:nrow(ora_res@result)]
  
  #we filter again since we fixed the FDR
  ora_res@result <-ora_res@result[ora_res@result$p.adjust <= pvalue,]
  #done
  return(ora_res)
}

#org --> "hsa" or "mmu"
KEGGenrichment <- function(
    df,
    org      = "mmu", 
    db       = org.Mm.eg.db, 
    pvalue   = 0.05, 
    qvalue   = 0.1,
    universe = NULL
) {
  require(clusterProfiler)
  library(KEGGREST)
  require(AnnotationDbi)
  # get pathways and their entrez gene ids
  
  path_entrez  <- keggLink("pathway", org) %>% 
    tibble(pathway = ., eg = sub(paste0(org,":"), "", names(.)))
  # get gene symbols and ensembl ids using entrez gene ids
  kegg_anno <- path_entrez %>%
    mutate(
      symbol = mapIds(db, eg, "SYMBOL", "ENTREZID"),
      ensembl = mapIds(db, eg, "ENSEMBL", "ENTREZID")
    )
  kegg_anno$pathway <- sub("path:", "", kegg_anno$pathway)
  # Pathway names
  pathways <- keggList("pathway", org) %>% 
    tibble(pathway = names(.), description = .)
  
  KEGG_pathways <- left_join(kegg_anno, pathways)
  KEGG_pathways$description <- sub(" - .*", "", KEGG_pathways$description )
  
  gs_kegg <- as.data.frame(KEGG_pathways[,c("description","symbol")])
  names(gs_kegg) <- c("term","gene")
  
  #now we use trick to fix clusterProfiler background problem https://github.com/YuLab-SMU/clusterProfiler/issues/283
  #trick from https://academic.oup.com/bioinformaticsadvances/article/4/1/vbae159/7829164#493480839
  # we add all genes to the universe as "background" term to avoid the problem
  if(!is.null(universe)){
    bgdf <- data.frame("background",universe)
    names(bgdf) <- c("term","gene")
    gs_kegg <- rbind(gs_kegg,bgdf)
  }  
  
  ora_res <- enricher(gene = df ,
                      universe = universe, 
                      minGSSize = 5, 
                      maxGSSize = 500000, 
                      TERM2GENE = gs_kegg,
                      pAdjustMethod="fdr", 
                      pvalueCutoff = 1,
                      qvalueCutoff = 1  )
  
  #add enrichment score
  gr <- as.numeric(sapply(strsplit(ora_res@result$GeneRatio,"/"),"[[",1)) /
    as.numeric(sapply(strsplit(ora_res@result$GeneRatio,"/"),"[[",2))
  
  br <- as.numeric(sapply(strsplit(ora_res@result$BgRatio,"/"),"[[",1)) /
    as.numeric(sapply(strsplit(ora_res@result$BgRatio,"/"),"[[",2))
  
  ora_res@result$es <- gr/br
  
  #now we fix the wrong FDR calculation
  nsets <- length(which(table(gs_kegg$term)>=5))
  nres <- nrow(ora_res@result)
  if (nres == 0) {
    return(data.frame())
  }
  diff <- nsets - nres
  pvals <- c(ora_res@result$pvalue,rep(1,diff))
  ora_res@result$p.adjust <- p.adjust(pvals,method="fdr")[1:nrow(ora_res@result)]
  #we filter again since we fixed the FDR
  ora_res@result <-ora_res@result[ora_res@result$p.adjust <= pvalue,]
  #done
  return(ora_res)
}


ReactomeEnrichment <- function(
    df,
    genome ,
    pvalue   = 0.05, 
    qvalue   = 0.1,
    universe = NULL
) {
  
  require(clusterProfiler)
  if(genome == "human"){
    pathToFile <- paste0('/shares/CIBIO-Storage/scratch_CIBIO/sharedLC/DavideBressan/RNA-seq/ALL_gene_sets/',genome,'/c2.cp.reactome.v2024.1.Hs.symbols.gmt')
  } else if(genome == "mouse"){
    pathToFile <- paste0('/shares/CIBIO-Storage/scratch_CIBIO/sharedLC/DavideBressan/RNA-seq/ALL_gene_sets/',genome,'/m2.cp.reactome.v2024.1.Mm.symbols.gmt')
  }
  
  stopifnot(file.exists(pathToFile))
  
  gs_general <- read.gmt(pathToFile)
  colnames(gs_general) <- c("term","gene")
  
  #now we use trick to fix clusterProfiler background problem https://github.com/YuLab-SMU/clusterProfiler/issues/283
  #trick from https://academic.oup.com/bioinformaticsadvances/article/4/1/vbae159/7829164#493480839
  # we add all genes to the universe as "background" term to avoid the problem
  if(!is.null(universe)){
    bgdf <- data.frame("background",universe)
    names(bgdf) <- c("term","gene")
    gs_general <- rbind(gs_general,bgdf)
  }  
  
  ora_res <- enricher(gene = df ,
                      universe = universe, 
                      minGSSize = 5, 
                      maxGSSize = 500000, 
                      TERM2GENE = gs_general,
                      pAdjustMethod="fdr", 
                      pvalueCutoff = 1,
                      qvalueCutoff = 1  )
  
  #add enrichment score
  gr <- as.numeric(sapply(strsplit(ora_res@result$GeneRatio,"/"),"[[",1)) /
    as.numeric(sapply(strsplit(ora_res@result$GeneRatio,"/"),"[[",2))
  
  br <- as.numeric(sapply(strsplit(ora_res@result$BgRatio,"/"),"[[",1)) /
    as.numeric(sapply(strsplit(ora_res@result$BgRatio,"/"),"[[",2))
  
  ora_res@result$es <- gr/br
  
  #now we fix the wrong FDR calculation
  nsets <- length(which(table(gs_general$term)>=5))
  nres <- nrow(ora_res@result)
  if (nres == 0) {
    return(data.frame())
  }
  diff <- nsets - nres
  pvals <- c(ora_res@result$pvalue,rep(1,diff))
  ora_res@result$p.adjust <- p.adjust(pvals,method="fdr")[1:nrow(ora_res@result)]
  #we filter again since we fixed the FDR
  ora_res@result <-ora_res@result[ora_res@result$p.adjust <= pvalue,]
  #done
  return(ora_res)
  
}

GeneralEnrichment <- function(
    df,
    genome ,
    pvalue   = 0.05, 
    qvalue   = 0.1,
    universe = NULL
) {
  
  require(clusterProfiler)
  if(genome == "human"){
    pathToFiles <- c(
      paste0('resources/',genome,'c3.tft.gtrd.v2024.1.Hs.symbols.gmt'),
      paste0('resources/',genome,'c8.all.v2024.1.Hs.symbols.gmt'),
      paste0('resources/',genome,'h.all.v2024.1.Hs.symbols.gmt'),
    )
  } else if(genome == "mouse"){
    pathToFiles <- c(
      paste0('resources/',genome,'/m3.gtrd.v2024.1.Mm.symbols.gmt'),
      paste0('resources/',genome,'/m8.all.v2024.1.Mm.symbols.gmt'),
      paste0('resources/',genome,'/mh.all.v2024.1.Mm.symbols.gmt'))
  }
  names(pathToFiles) <- c("TF targets - GTRD","Cell Types","Hallmark")
  
  lapply(pathToFiles, function(x) stopifnot(file.exists(x)))
  
  allSets <- lapply(pathToFiles, function(x){
    tmpDf <- read.gmt(x)
    colnames(tmpDf) <- c("term","gene")
    tmpDf
  })
  
  
  allResults <- lapply(allSets, function(gs_general) {
    #now we use trick to fix clusterProfiler background problem https://github.com/YuLab-SMU/clusterProfiler/issues/283
    #trick from https://academic.oup.com/bioinformaticsadvances/article/4/1/vbae159/7829164#493480839
    # we add all genes to the universe as "background" term to avoid the problem
    if(!is.null(universe)){
      bgdf <- data.frame("background",universe)
      names(bgdf) <- c("term","gene")
      gs_general <- rbind(gs_general,bgdf)
    }  
    
    ora_res <- enricher(gene = df ,
                        universe = universe, 
                        minGSSize = 5, 
                        maxGSSize = 500000, 
                        TERM2GENE = gs_general,
                        pAdjustMethod="fdr", 
                        pvalueCutoff = 1,
                        qvalueCutoff = 1  )
    
    #add enrichment score
    gr <- as.numeric(sapply(strsplit(ora_res@result$GeneRatio,"/"),"[[",1)) /
      as.numeric(sapply(strsplit(ora_res@result$GeneRatio,"/"),"[[",2))
    
    br <- as.numeric(sapply(strsplit(ora_res@result$BgRatio,"/"),"[[",1)) /
      as.numeric(sapply(strsplit(ora_res@result$BgRatio,"/"),"[[",2))
    
    ora_res@result$es <- gr/br
    #now we fix the wrong FDR calculation
    nsets <- length(which(table(gs_general$term)>=5))
    nres <- nrow(ora_res@result)
    if (nres == 0) {
      return(data.frame())
    }
    diff <- nsets - nres
    pvals <- c(ora_res@result$pvalue,rep(1,diff))
    ora_res@result$p.adjust <- p.adjust(pvals,method="fdr")[1:nrow(ora_res@result)]
    #we filter again since we fixed the FDR
    ora_res@result <-ora_res@result[ora_res@result$p.adjust <= pvalue,]
    #done
    return(ora_res)
  })
  return(allResults)
}



GSEA_enrichment <- function(df,genome, pathways.gmt) {
  if(genome == "human"){
    if(pathways.gmt == "Hallmark"){
      pathways <- paste0('resources/',genome,'/h.all.v2024.1.Hs.symbols.gmt')
    } else if(pathways.gmt == "C2/M2"){
      pathways <- paste0('resources/',genome,'/c2.all.v2024.1.Hs.symbols.gmt')
    } else if(pathways.gmt == "C3/M3"){
      pathways <- paste0('resources/',genome,'/c3.tft.gtrd.v2024.1.Hs.symbols.gmt')
    }
  } else if(genome == "mouse"){
    if(pathways.gmt == "Hallmark"){
      pathways <- paste0('resources/',genome,'/mh.all.v2024.1.Mm.symbols.gmt')
    } else if(pathways.gmt == "C2/M2"){
      pathways <- paste0('resources/',genome,'/m2.all.v2024.1.Mm.symbols.gmt')
    } else if(pathways.gmt == "C3/M3"){
      pathways <- paste0('resources/',genome,'/m3.gtrd.v2024.1.Mm.symbols.gmt')
    }
  }
  pathwaysSet <- fgsea::gmtPathways(pathways)
  
  # Create a list containing a named vector (with genenames) of log2fc
  geneList <- df$log2FoldChange
  names(geneList) <- toupper(df$Geneid)
  geneList <- geneList[!is.na(geneList)]
  # Run GSEA algorithm
  fgseaRes <- fgsea::fgsea(pathways = pathwaysSet, 
                           stats   = geneList,
                           minSize = 10,
                           nperm = 1000,
                           maxSize = 2000)
  if (nrow(fgseaRes) == 0) {
    return(data.frame())
  } else {
    return(fgseaRes)
  }
}


#----------------------------------------------------------------------------------------------
# Volcano plot code
#----------------------------------------------------------------------------------------------
VolcanoPlot <- function(df, xlim=NULL, ylim=NULL, main = NULL, labelSize = 8, pval = 0.05, log2FC = 1) {
  require(ggplot2)
  require(dplyr)
  # require(ggrastr)

  df <- mutate(df, shape = "circle")

  df <- mutate(df, shape = ifelse(-log10(padj) >  ylim[2], "triangle", shape))
  df <- mutate(df, padj = ifelse(-log10(padj) >  ylim[2], 10^-ylim[2], padj))

  df <- mutate(df, shape = ifelse(log2FoldChange > xlim[2], "triangle", shape))
  df <- mutate(df, shape = ifelse(log2FoldChange < -xlim[2], "triangle", shape))

  df <- mutate(df, log2FoldChange = ifelse(log2FoldChange > xlim[2], xlim[2], log2FoldChange))
  df <- mutate(df, log2FoldChange = ifelse(log2FoldChange < -xlim[2], -xlim[2], log2FoldChange))


  p <-  ggplot(data = na.omit(df), aes(x=log2FoldChange, y=-log10(padj), colour=DEG, shape=shape) ) +

    # geom_point_rast(alpha=0.7, size=1.7, raster.height = 5.15, raster.width = 6, raster.dpi = 400) +
    geom_point(alpha=0.7, size=1.7) +

    annotate("text", label = sum(df$DEG == "Upregulated"), color = "red", y = 0, x = xlim[2],
             vjust="inward",hjust="inward", size = labelSize) +
    annotate("text", label = sum(df$DEG == "Downregulated"), color = "darkgreen", y = 0, x = xlim[1],
             vjust="inward",hjust="inward", size = labelSize) +

    theme_classic(base_size = 20) +
    theme(legend.title = element_blank()) +
    theme(legend.position = "top") +

    ggtitle(main) +
    theme(plot.title = element_text(lineheight=.8, face="bold", hjust = .5)) +

    xlim(xlim) + ylim(ylim) +

    geom_hline(yintercept = -log10(pval), linetype = 2) +
    geom_vline(xintercept = c(-log2FC, log2FC), linetype = 2) +

    xlab(bquote(~Log[2]~ "fold change")) + ylab(bquote(~-Log[10]~italic(P))) +

    scale_colour_manual(values=c("Downregulated" = "darkgreen", "NS" = "gray", "Upregulated" = "red"),
                        labels = c("Downregulated" = "Downregulated", "NS" = "NS", "Upregulated" = "Upregulated"),
                        drop = FALSE) + #Force legend to show always

    guides(shape=FALSE) # Remove legend for shapes

  return(p)
}



# VolcanoPlot <- function(df, xlim=NULL, ylim=NULL, main = NULL, labelSize = 8, pval = 0.05, log2FC = 1) {
#   require(ggplot2)
#   require(ggrastr)

#   df <- mutate(df, shape = "circle")
#   df <- mutate(df, shape = ifelse(-log10(padj) >  ylim[2], "triangle_up", shape))
#   df <- mutate(df, padj = ifelse(-log10(padj) >  ylim[2], 10^-ylim[2], padj))

#   df <- mutate(df, shape = ifelse(log2FoldChange > xlim[2], "triangle_right", shape))
#   df <- mutate(df, shape = ifelse(log2FoldChange < -xlim[2], "triangle_left", shape))

#   df <- mutate(df, log2FoldChange = ifelse(log2FoldChange > xlim[2], xlim[2], log2FoldChange))
#   df <- mutate(df, log2FoldChange = ifelse(log2FoldChange < -xlim[2], -xlim[2], log2FoldChange))


#   p <-  ggplot(data = na.omit(df), aes(x=log2FoldChange, y=-log10(padj), colour=DEG, shape=shape) ) +

#     # geom_point_rast(alpha=0.7, size=1.7, raster.height = 5.15, raster.width = 6, raster.dpi = 400) +
#     geom_point(alpha=0.7, size=1.7) +
    
#     annotate("text", label = sum(df$DEG == "Upregulated"), color = "red", y = 0, x = xlim[2],
#              vjust="inward",hjust="inward", size = labelSize) +
#     annotate("text", label = sum(df$DEG == "Downregulated"), color = "darkgreen", y = 0, x = xlim[1],
#              vjust="inward",hjust="inward", size = labelSize) +

#     theme_classic(base_size = 20) +
#     theme(legend.title = element_blank()) +
#     theme(legend.position = "top") +


#     ggtitle(main) +
#     theme(plot.title = element_text(lineheight=.8, face="bold", hjust = .5)) +

#     xlim(xlim) + ylim(ylim) +

#     geom_hline(yintercept = -log10(pval), linetype = 2) +
#     geom_vline(xintercept = c(-log2FC, log2FC), linetype = 2) +

#     xlab(bquote(~Log[2]~ "fold change")) + ylab(bquote(~-Log[10]~italic(P))) +

#     scale_colour_manual(values=c("Downregulated" = "darkgreen", "NotDE" = "gray", "Upregulated" = "red"),
#                         labels = c("Downregulated" = "Downregulated", "NS" = "NS", "Upregulated" = "Upregulated"),
#                         drop = FALSE) + #Force legend to show always

#     scale_shape_manual(values=c("triangle_up" = "\u25B2", "triangle_right" = "\u25BA", "triangle_left" = "\u25C4", "circle" = "\u25CF",
#                                 drop = FALSE)) +

#     guides(shape=FALSE) # Remove legend for shapes

#   return(p)
# }