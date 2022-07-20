runtopGO <- function(topGOdata, topNodes = 20, conservative = TRUE) {
  
  RFisher <- runTest(topGOdata, 
    algorithm = "classic", 
    statistic = "fisher", verbose = F)
  
  # To make this test conservative. Next we will test the enrichment using the Kolmogorov-Smirnov test. We will use the both the classic and the elim method.
  
  if(conservative) 
  {
    RKS <- runTest(topGOdata, algorithm = "classic", 
      statistic = "ks", verbose = F)
    
    RKS.elim <- runTest(topGOdata, algorithm = "elim", 
      statistic = "ks", verbose = F)
    
    
    
    
    allRes <- GenTable(topGOdata, 
      classicFisher = RFisher,
      classicKS = RKS, 
      elimKS = RKS.elim,
      orderBy = "elimKS", 
      ranksOf = "classicFisher", 
      topNodes = topNodes) 
  } else {
    RKS <- runTest(topGOdata, algorithm = "classic", 
      statistic = "ks")
    
    test.stat <- new("weightCount",
      testStatistic = GOFisherTest,
      name = "Fisher test", sigRatio = "ratio")
    
    weights <- getSigGroups(topGOdata, test.stat)
    
    allRes <- GenTable(topGOdata,
      classic = RFisher,
      KS = RKS,
      weight = weights,
      orderBy = "weight",
      ranksOf = "classic",
      topNodes = topNodes)
    
    # allRes <- GenTable(object = topGOdata, 
    #                    elimFisher = RFisher,
    #                    topNodes = topNodes)
  }
  
  return(allRes)
}

GOenrichment <- function(df, cons = T, onto = "BP", Nodes = Inf) {
  
  # df needs to include follow cols: "ids"            "sampleB"        "GO.ID" (list)      and    "moduleMemberCorr"
  # 
  # moduleMemberCorr
  
  df %>% pull(moduleMemberCorr) -> query.p 
  df %>% pull(ids) -> query.names
  
  names(query.p) <- query.names
  
  
  
  # keep MAP of query genes
  # gene2GO <- df$GO.ID
  # names(gene2GO) <- query.names
  
  # gene2GO <- MAP # unlist(MAP)
  gene2GO <- MAP[names(MAP) %in% query.names]
  
  cat("\n")
  cat("Annotation : ",length(gene2GO))
  cat("\n")
  cat("GO terms : ", length(unlist(gene2GO)))
  cat("\n")
  
  description <- "complete topGO enrichment using split_annot"
  
  topCor <- function(x) { return(abs(x) > 0.5) }
  
  topGOdata <- new("topGOdata", 
    ontology = onto, 
    description = description,
    allGenes = query.p,
    # geneSel: function to specify which genes are interesting based on the gene scores. 
    # It should be present if the allGenes object is of type numeric.
    # geneSel = function(x) { x == 1 }, # for log-fc
    geneSel = topCor,
    # geneSel = function(x) x, # for p value
    annot = annFUN.gene2GO,
    # mapping = hsGO, # omit this flag
    gene2GO = gene2GO)
  
  # run TopGO results 
  
  allGO <- usedGO(topGOdata)
  
  if(is.infinite(Nodes)) {
    topNodes <- length(allGO)
  } else {
    topNodes <- Nodes
  }
  
  
  allRes <- runtopGO(topGOdata, topNodes = Nodes, conservative = cons)
  
  # make p adjustable
  
  p.adj.ks <- p.adjust(allRes$classicKS , method="BH")
  
  allRes <- cbind(allRes, p.adj.ks)
  
  allRes$Term <- gsub(" [a-z]*\\.\\.\\.$", "", allRes$Term)
  allRes$Term <- gsub("\\.\\.\\.$", "", allRes$Term)
  
  return(allRes)
  
  
}

prep_dist_data <- function(df, select_exclusive = T) {
  
  # names(df) should contain at least cols:
  # "ids"            
  # "sampleB"       
  #  "log2FoldChange" 
  # "pvalue" or    "padj" 
  # "lfcT"
  
  
  if(select_exclusive) {
    
    # Select only n_sam == 1 ----
    
    df %>% # this is the only input needed
      dplyr::select(ids, sampleB) %>% 
      group_by(ids) %>%
      summarise(n_sam = length(sampleB), 
        across(sampleB, .fns = list)) -> distinct_trans
    
    # Sanity check of correct grouping exclusive ids
    # Also you can compare against the upset png value bars
    
    distinct_trans %>% group_by(n_sam) %>% tally() -> prevalence
    
    distinct_trans %>% 
      filter(n_sam == 1) %>% 
      mutate(sampleB = unlist(sampleB)) -> distinct_trans
    
    cat("\n Prevalence: ",
      paste0(prevalence$n_sam, 's'), "\n", 
      "\t",prevalence$n, "\n")
    
  } else {
    df %>% dplyr::select(ids, sampleB) %>% 
      group_by(sampleB) %>% distinct(ids) -> distinct_trans
    
    distinct_trans %>% group_by(sampleB) %>% tally() -> prevalence
    
    cat("\n Prevalence: ",
      paste0(prevalence$sampleB, 's'), "\n", 
      "\t",prevalence$n, "\n")
    
  }
  
  
  
  
  
  # 2) Subset of GO ids based on the distinct transcript dataset ----
  
  distinct_trans %>% pull(ids) %>% unique(.) -> query.ids
  
  n <- sum(keep <- names(MAP) %in% query.ids) / length(query.ids) 
  
  cat('\n % of ids mapped in Gen Ontology: ', n*100, '%\n')
  
  STRG2GO <- data.frame(ids = rep(names(MAP[keep]),
    sapply(MAP[keep], length)),
    GO.ID = unlist(MAP[keep]), row.names = NULL) %>% as_tibble()
  
  # Sanity check: % of recovery querys with go ids
  # sum(unique(STRG2GO$ids) %in% query.ids) / length(query.ids) 
  
  distinct_trans %>% 
    left_join(STRG2GO, by = "ids") %>%
    group_by(ids, sampleB) %>%
    summarise(across(GO.ID, .fns = list), .groups = "drop_last") %>%
    ungroup() -> distinct_trans
  
  # Add p values to the distinct trans
  distinct_trans %>% 
    left_join(df %>% dplyr::select(ids, padj, log2FoldChange), 
      by = "ids") -> distinct_trans
  
  
  return(distinct_trans)
  
}

# anntation

split_blast <- function (x, hit = "sprot_Top_BLASTX_hit")
{
  require(data.table)
  y <- x[!is.na(get(hit)), .(get(hit), gene_id, transcript_id,
    prot_id)]
  z <- strsplit(y$V1, "`")
  n <- sapply(z, length)
  z <- strsplit(unlist(z), "\\^")
  if (any(sapply(z, "[", 1) != sapply(z, "[", 2)))
    print("WARNING: check different values in columns 1 and 2")
  NAME <- gsub("^RecName: Full=", "", sapply(z, "[", 6))
  NAME <- gsub("SubName: Full=", "", NAME)
  NAME <- gsub(";$", "", NAME)
  NAME <- gsub(" \\{[^}]+}", "", NAME)
  x1 <- data.frame(gene = rep(y$gene_id, n), transcript = rep(y$transcript_id,
    n), protein = rep(gsub(".*\\|", "", y$prot_id), n), uniprot = sapply(z,
      "[", 1), align = sapply(z, "[", 3), identity = as.numeric(gsub("%ID",
        "", sapply(z, "[", 4))), evalue = as.numeric(gsub("E:",
          "", sapply(z, "[", 5))), name = NAME, lineage = sapply(z,
            "[", 7), domain = gsub("; .*", "", sapply(z, "[", 7)),
    genus = gsub(".*; ", "", sapply(z, "[", 7)), stringsAsFactors = FALSE)
  message(nrow(x1), " ", hit, " annotations")
  data.table(x1)
}

# functions to test, check profiling_oyster_set.R script

# enrichement by GO/entrezid

getEntrez <- function(go, de_df, pattern = NULL, GOtype='MF', mart) {
  
  library(biomaRt)
  
  if(is.null(pattern)) {
    geneList <- getList(go, de_df, '', GOtype)
  } else
    geneList <- getList(go, de_df, pattern, GOtype)
  
  
  
  GO_universe <- geneList %>% select_at(vars(go)) %>% pull(.)
  
  goids = getBM(attributes = c('entrezgene_id', 
    'go_id',
    'name_1006'),
    filters = 'go', 
    values = GO_universe, 
    mart = mart)
  
  
  entrezgene_id <- goids %>%
    filter(!is.na(entrezgene_id)) %>%
    mutate(entrezgene_id = as.character(entrezgene_id)) %>%
    arrange(desc(entrezgene_id))
  
  return(entrezgene_id)
}

# Example
# getEntrez(go, de_df, GOtype='MF', mart = ensembl_cgigas,
#   pattern = NULL) -> mart_data
# 

# wgcna data viz

plotscatterSig <- function(x, M = ..., sig = 0.05) {
  
  colVal <- c("#54278f", "grey30")
  
  x %>%
    mutate(color = ifelse(p.GS < sig, '*', 'ns')) %>%
    filter(module %in% M) %>%
    # filter(group != 'HC') %>%
    ggplot(aes(x = moduleMemberCorr, y = abs(GS), color = color, alpha = -log10(p.GS))) +  
    # geom_vline(color = "black", xintercept = 0.8, linetype="dashed", alpha=0.5) +
    # geom_hline(color = "black", yintercept = 0.8, linetype="dashed", alpha=0.5) +
    # geom_hline(color = "black", yintercept = -0.8, linetype="dashed", alpha=0.5) +
    geom_point(size = 1) +
    labs(x = 'Module membership', y = 'Gene significance') + # 'Gene Correlation'
    scale_color_manual(name = '', values = colVal) + # breaks = c(0, 0.25 ,0.5, 0.75, 1), limits= c(0, 1)
    scale_alpha_continuous(name = expression(-Log[10] ~ "P")) +
    facet_grid(module~group) +
    # ggh4x::facet_nested_wrap(module+group ~., nest_line = T) +
    theme_bw(base_family = "GillSans", base_size = 16) +
    # geom_smooth(method = "lm", linetype="dashed", size = 0.5, alpha=0.5, color = 'white',
    #   se = TRUE, na.rm = TRUE) +
    theme(legend.position = 'top') -> psave
  
  filename <- paste0('geneTraitSignificance_', paste(M, collapse = '_'), '.png')
  
  ggsave(psave, filename = filename, path = path, width = 10, height = 5) 
  
  # psave
  
}

# egg


get_eggnog <- function (x, ids, by = "transcript_id") 
{
  # trinotateR::plot_NOGs()
  
  names(x)[grepl('transcript', names(x))] <- 'transcript_id'
  names(x)[grepl('egg', names(x))] <- 'eggnog'
  
  if (by == "transcript_id") {
    
    x1 <- x %>% filter(transcript_id %in% ids)
    
    # y <- unique(x1[!is.na(eggnog), .(transcript_id, eggnog)])
    eggnog <- x1 %>% drop_na(eggnog) %>% distinct(eggnog) %>% pull(eggnog)
    
  }
  else {
    x1 <- x %>% filter(gene_id %in% ids)
    
    y <- unique(x1[!is.na(eggnog), .(gene_id, eggnog)])
    
  }
  # nogs <- gsub("(.*)\\^.*", "\\1", eggnog)
  
  nogs <- eggnog
  
  
  # y %>% separate(col = eggnog, sep = "(.*)\\^.*", into = c('nogs', 'eggnog'))
  
  f <- list.files(pattern = 'NOG.annotations.tsv', path = getwd())
  
  if(identical(f, character(0))) {
    
    download.file("https://raw.githubusercontent.com/RJEGR/infovis/master/NOG.annotations.tsv", "NOG.annotations.tsv")
    
    egg <- read.table("NOG.annotations.tsv", sep="\t", stringsAsFactors=FALSE, quote="")
    
  } else {
    egg <- read.table(f , sep="\t", stringsAsFactors=FALSE, quote="")
  }
  
  names(egg) <- c("db", "nog", "proteins", "species", "class", "description")
  
  n <- match(nogs, egg$nog)
  
  y <- table(unlist(strsplit(egg$class[n], "")))
  
  y <- data.frame(y)
  
  names(y) <- c('code', 'Freq')
  
  
  
  return(y)
}

split_KEGG <- function (x, hit = "Kegg")  {
  library(data.table)
  y <- x[!is.na(get(hit)), .(get(hit), gene_id, transcript_id, 
    prot_id)]
  z <- strsplit(y$V1, "`")
  n <- sapply(z, length)
  z <- strsplit(unlist(z), "\\^")
  x1 <- data.frame(gene = rep(y$gene_id, n), 
    transcript = rep(y$transcript_id, n), 
    protein = rep(gsub(".*\\|", "", y$prot_id), n), 
    Kegg = sapply(z, "[", 1), stringsAsFactors = FALSE)
  message(nrow(x1), " ", hit, " annotations")
  data.table(x1)
}


