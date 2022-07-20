
library(ggraph)
library(igraph)
library(tidygraph)
library(WGCNA)
library(tidyverse)


rm(list = ls()) # Limpiar la memoria de la sesion de R
if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE) # 


path <- "~/Downloads/CORAL/"

egg_df <- read_rds(paste0(path, 'trinotate_eggnog.rds'))
blastp <- readRDS(paste0(path, 'blastp.rds'))
allData_df <- readRDS(paste0(path,'allData_df.rds'))
load(paste0(path, 'moduleTraitRelationship.RData'))

# Select query transcripts based on significance modules 
allData_df %>% filter(module %in% psME) %>% pull(ids) %>% unique() -> query.transcripts


sum(allData_df$ids %in% query.transcripts) / length(allData_df$ids) # 60 % of total transcriptome is kept


# (Omit) ----
# this step is only one time to create Nodes and edges txt files 

bwnet <- readRDS(paste0(path, "bwnet.rds"))

bwModuleColors = labels2colors(bwnet$colors)

names(bwModuleColors) <- names(bwnet$colors)

data.frame(module = bwModuleColors, 
  blocks = bwnet$blocks) %>% 
  as_tibble(rownames = 'transcripts') %>%
  filter(transcripts %in% query.transcripts) %>%
  group_by(module, blocks) %>%
  dplyr::count() %>%
  pull(blocks) -> blocksN


exportNet <- function(TOMobj, blockGenes, moduleColors, threshold = 0.9) {
  
  nodeNames <- names(moduleColors)[blockGenes]
  nodeAttr <- moduleColors[blockGenes]
  
  
  file_out <- gsub('.RData', '', basename(TOMFile))
  
  edgeFile <- paste0(path, "/",file_out, ".edges.txt")
  nodeFile <- paste0(path, "/",file_out, ".nodes.txt")
  
  
  Net = exportNetworkToCytoscape(TOMobj,
    edgeFile = edgeFile,
    nodeFile = nodeFile,
    weighted = TRUE,
    threshold = threshold,
    nodeNames = nodeNames, nodeAttr = nodeAttr)
  
  cat('\nEdges: ',nrow(Net$edgeData), 'and Nodes: ',nrow(Net$nodeData),'\n')
  # dim(Edges <- Net$edgeData)
  # dim(Nodes <- Net$nodeData)
  
  graph <- tbl_graph(nodes = Net$nodeData, 
    edges = Net$edgeData, directed = FALSE)
  
  return(graph)
  
  # if(is.null(query)) {
  #   graph
  # } else {
  #   graph %>% activate(nodes) %>% filter(nodeName %in% query)
  # }
  
  
  # return(Net)
  
}

Nets <- list()

n <- length(bwnet$TOMFiles)

for(i in n:1) {
  j <- i
  
  TOMFile <- list.files(path = path, pattern = bwnet$TOMFiles[j], full.names = T)
  
  cat('\nReading', basename(TOMFile), '\n')
  
  blockGenes <- bwnet$blockGenes[[j]]
  
  load(TOMFile)
  
  TOMobj <- TOM
  
  exportNet(TOMobj, blockGenes, bwModuleColors, threshold = 0.7) -> Net
  
  Nets[[j]] <- Net
}


# Continue ----

EdgeFiles <- list.files(path = path, pattern = ".edges.txt", full.names = T)
NodesFiles <- list.files(path = path, pattern = ".nodes.txt", full.names = T)

Edges <- do.call(rbind, lapply(EdgeFiles, data.table::fread))
Nodes <- do.call(rbind, lapply(NodesFiles, data.table::fread))


names(Nodes)[3] <- 'module'

# table(Nodes$module)

Nodes <- Nodes %>% select(-altName)
# Sanity check
length(Nodes$nodeName) # 27252 transcripts w/ threshold > 0.7 during exportNet step
length(query.transcripts) # 28950

sum(Nodes$nodeName %in% query.transcripts) # 23094 istead of 28950 !?

# Annotation features ----

# Enrich first nodes by domain and genus name from blastp (unique names) 
# Then, enrich by the GO.ID (Gene Ontology) 

# By blastp ----

length(unique(Nodes$nodeName) -> query.names)

blastp <- readRDS(paste0(path, 'blastp.rds')) %>% 
  select(transcript, domain, genus, identity, evalue) %>% 
  filter(transcript %in% query.names) %>%
  rename('nodeName' = 'transcript')


# blastp %>% count(domain, sort = T) 

# blastp %>% 
#   filter(domain %in% 'Eukaryota') %>%
#   # mutate(genus = fct_reorder(genus, identity)) %>%
#   count(genus, sort = T)

# blastp %>% 
#   filter(domain %in% 'Eukaryota') %>%
#   mutate(genus = fct_reorder(genus, identity)) %>%
#   # filter(identity > 50) %>%
#   ggplot(aes(y = identity, x = genus)) +
#   geom_boxplot() +
#   coord_flip()

nodeInfo <- blastp %>% 
  distinct(nodeName, genus, .keep_all = T) %>% 
  select(-identity,-evalue) %>%
  right_join(Nodes) %>%
  mutate(genus = ifelse(is.na(genus), '', genus),
    domain = ifelse(is.na(domain), '', domain))

# Gene Ontology  ---- 
# Formular el metodo correcto para analizar nodos en base a GO, pues se tiene mas de un go.id para muchos nodos

llist <- function(x) {
  x <- paste(x, sep = ';', collapse = ';')
  x <- list(x)
  # x <- unlist(x)
}

.go <- readRDS(paste0(path, 'gene_ontology_BLASTP.rds')) %>%
  filter(ontology %in% 'biological_process') %>%
  filter(transcript %in% query.names) %>%
  rename('nodeName' = 'transcript')

go <- .go %>%
  # select(-) %>%
  group_by(nodeName) %>%
  summarise(go_freq = length(go), 
    across(go, .fns = llist))

nodeInfo <- go %>% right_join(nodeInfo) 

# Test term sem

nodeInfo %>% 
  filter(go_freq == 1) %>%
  mutate(go = unlist(go)) %>%
  left_join(.go)  -> .nodeInfo

library(rrvgo)

sem <- function(x) {
  # if(n > 1) {
  hsGO <- read_rds(paste0(path, '/hsGO_BP.rds'))
  
  # x <- nodeInfo[1,3]$go
  x <- strsplit(unlist(x), ";")[[1]]
  
  SimMatrix <- calculateSimMatrix(x, 
    orgdb = 'org.Hs.eg.db',
    ont="BP", 
    semdata = hsGO,
    method = 'Wang')
  
  reduceSimMatrix(SimMatrix, threshold = 0.9, orgdb = 'org.Hs.eg.db')
  # } else
  
}

# Works good!!

nodeInfo %>% 
  filter(go_freq != 1) %>%
  head() %>%
  group_by(nodeName) %>%
  summarise(sem(go))

# go_file <- paste0(path, 'Trinotate_report.xls.gene_ontology')
# 
# MAP <- topGO::readMappings(go_file)

# Prepare igraph ----

graph <- tbl_graph(nodes = nodeInfo, edges = Edges, directed = FALSE)

saveRDS(graph, file = paste0(path, 'graph.rds'))

str(query.transcripts[query.transcripts %in% nodeInfo$nodeName] -> query.transcripts)

# Filtering nodes and edges

graph %>%  activate("nodes") %>% filter(module %in% psME) -> g

write_rds(g, file = paste0(path, paste(psME, collapse = '_'), '_graph.rds'))
