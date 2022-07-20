
library(ggraph)
library(igraph)
library(tidygraph)
library(tidyverse)


rm(list = ls()) # Limpiar la memoria de la sesion de R
if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE) # 

psME <- c('blue', 'turquoise')

path <- "~/Downloads/CORAL/"

fileName <- "blue_turquoise_graph.rds"

g <- read_rds(file = paste0(path, fileName))

# Sanity check

g %>% activate("nodes") %>% as_tibble() %>% count(module)

# blue       9144
# turquoise 13976

# if table(Nodes$module)[names(table(Nodes$module)) %in% psME] 

# blue 9133
# turquoise 13961 

g %>% activate("nodes") %>%  
  mutate(hub = centrality_hub(),
    degree = centrality_degree(),
    eigen = centrality_eigen(),
    pageR = centrality_pagerank()) -> g

col_palette <- structure(psME, names = psME)

library(rstatix)

g %>% activate("nodes") %>% as_tibble() %>% 
  cor_mat(vars = c('hub', 'degree', 'eigen', 'pageR'), method = 'spearman') %>%
  cor_reorder() %>%
  pull_lower_triangle() %>% 
  cor_plot(label = TRUE, method = 'color', insignificant = "blank")
  # cor_get_pval()

# Hub and eigen as a highest relation to centrality degree

# Further reads https://cambridge-intelligence.com/eigencentrality-pagerank/

g %>% activate("nodes") %>% as_tibble() %>% 
  ggplot(aes(degree, color = module)) +
  # facet_grid(module ~ name, scales = 'free') +
  geom_freqpoly(size = 1, binwidth = 50) +
  theme_bw(base_family = "GillSans", base_size = 14) +
  scale_color_manual('Module', values = col_palette) +
  labs(y = 'Transcripts', x = 'Centrality degree') +
  theme(legend.position = 'top',
    strip.background = element_rect(fill = 'grey86', color = 'white'),
    panel.border = element_blank()) -> psave

filename <- paste0('degree_', paste(psME, collapse = '_'), '.png')


ggsave(psave, path = path, filename = filename, 
  width = 5, height = 5, dpi = 1000) #

degree_f <- function(x) {quantile(x, probs = 0.95)}

g %>% activate('nodes') %>% 
  group_by(module) %>%
  mutate(nodeName = ifelse(degree > degree_f(degree), nodeName, NA)) -> g
# filter(degree > quantile(V(.)$degree, probs = 0.95)) -> g

g %>% activate("nodes") %>% as_tibble() %>% drop_na(nodeName) %>% count(module)

# blue        458
# turquoise   696

g %>% activate("nodes") %>% 
  ungroup() %>%
  filter(!is.na(nodeName)) -> .g

.g %>% activate("nodes") %>% as_tibble() %>% 
  cor_mat(vars = c('hub', 'degree', 'eigen', 'pageR'), method = 'spearman') %>%
  cor_reorder() %>%
  pull_lower_triangle() 

# g %>% ungroup() %>%
#   filter(!is.na(nodeName)) %>%
#   mutate(hub = centrality_hub(),
#       degree = centrality_degree(),
#       eigen = centrality_eigen(),
#   pageR = centrality_pagerank()) -> .g

hist(E(.g)$weight)

.g %>% activate("nodes") %>% as_tibble() %>% 
  pivot_longer(cols = c('hub', 'eigen', 'pageR')) %>%
  # filter(value > 0) %>%
  # filter(module %in% 'blue') %>%
  ggplot(aes(degree, value, color = module)) +
  # facet_grid(module ~ name, scales = 'free') +
  ggh4x::facet_nested_wrap(~module +  name, scales = 'free') +
  geom_jitter(alpha = 0.5, shape = 1) +
  scale_color_manual('Module', values = col_palette) +
  theme_bw(base_family = "GillSans", base_size = 14) +
  # labs(y = 'Transcripts', x = 'Centrality degree') +
  theme(legend.position = 'top',
    strip.background = element_blank(),
    panel.border = element_blank(),
    ggh4x.facet.nestline = element_line(colour = "black")) -> p

ggsave(p, path = path, filename = 'transcripts_centrality.png', 
  width = 12, height = 7, dpi = 1000) #


w_filter <- quantile(E(.g)$weight, probs = c(0.9))

.g %>% activate("edges") %>%
  mutate(weight = ifelse(weight > w_filter, weight, NA)) -> g

g %>% activate("edges") %>% filter(!is.na(weight))


hist(E(g)$weight)

# g %>% activate("nodes") %>% 
#   ungroup() %>%
#   filter(!is.na(nodeName)) -> g
# 
# .g %>% activate("edges") %>% 
#   mutate(edge_betweenness = centrality_edge_betweenness(weights = weight, directed = FALSE)) %>%
#   arrange(desc(edge_betweenness)) -> .g

# summary(E(.g)$edge_betweenness)

# hist(E(.g)$edge_betweenness)

# edge_filter <- quantile(E(g)$edge_betweenness, probs = 0.95)

# g %>% activate("edges") %>%
#   filter(edge_betweenness > edge_filter) -> g

# 
# g %>% activate("nodes") %>%  
#   # mutate(degree = centrality_degree()) %>%
#   filter(degree > 0) # -> g
# 


# continue ----
# 

g %>% activate("nodes") %>%  
  mutate(
    louvain = igraph::cluster_louvain(.)$membership,
    walktrap = cluster_walktrap(.)$membership,
    fast_greedy = cluster_fast_greedy(.)) -> g

hist(V(g)$louvain)
hist(V(g)$walktrap)
hist(V(g)$fast_greedy)

table(V(g)$domain)

g %>% activate("nodes") %>% as_tibble() -> nodeInfo



# go enrichment of nodes ----
nodeInfo %>% pull(nodeName) -> query.transcripts

llist <- function(x) {
  x <- paste(x, sep = ';', collapse = ';')
  x <- list(x)
  # x <- unlist(x)
}

.go <- readRDS(paste0(path, 'gene_ontology_BLASTP.rds')) %>%
  filter(ontology %in% 'biological_process') %>%
  filter(transcript %in% query.transcripts) %>%
  rename('nodeName' = 'transcript')

go <- .go %>%
  # select(-) %>%
  group_by(nodeName) %>%
  summarise(go_freq = length(go), 
    across(go, .fns = llist))

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
  if(is.list(x)) {
    
    x <- strsplit(unlist(x), ";")[[1]]
    
    x <- sort(x)
    
    SimMatrix <- calculateSimMatrix(x, 
      orgdb = 'org.Hs.eg.db',
      ont="BP", 
      semdata = hsGO,
      method = 'Wang')
    
    data <- reduceSimMatrix(SimMatrix, threshold = 0.9, orgdb = 'org.Hs.eg.db')

    data %>%
      arrange(desc(size)) %>%
      summarise(parent_freq = length(unique(parent)),
        parentTerm = head(.,1)$parentTerm,
        size = head(.,1)$size,
        across(parent, .fns = llist)) -> data

    
  } else {
    
    x <- sort(x)
    
    SimMatrix <- calculateSimMatrix(x, 
      orgdb = 'org.Hs.eg.db',
      ont="BP", 
      semdata = hsGO,
      method = 'Wang')
    
    reducedTerms <- reduceSimMatrix(SimMatrix, threshold = 0.9, orgdb = 'org.Hs.eg.db')
    
    
    y <- cmdscale(as.matrix(as.dist(1 - SimMatrix)), eig = TRUE, k = 2)
    
    data <- cbind(as.data.frame(y$points), 
      reducedTerms[match(rownames(y$points), reducedTerms$go),])
    
    # return(data)
  }

  # 
  # SimMatrix <- calculateSimMatrix(x, 
  #   orgdb = 'org.Hs.eg.db',
  #   ont="BP", 
  #   semdata = hsGO,
  #   method = 'Wang')
  # 
  # reducedTerms <- reduceSimMatrix(SimMatrix, threshold = 0.9, orgdb = 'org.Hs.eg.db')
  # 
  # 
  # y <- cmdscale(as.matrix(as.dist(1 - SimMatrix)), eig = TRUE, k = 2)
  # 
  # data <- cbind(as.data.frame(y$points), 
  #   reducedTerms[match(rownames(y$points), reducedTerms$go),])
  # 
  return(data)
  
}


nodeInfo %>%
  filter(go_freq != 1) %>%
  arrange(desc(go_freq)) %>%
  # head(1) %>%
  group_by(nodeName) %>%
  summarise(sem(go)) -> parentTerm

# Works good!!

# Single nodes ----

nodeInfo %>% 
  filter(go_freq == 1) %>%
  mutate(go = unlist(go)) %>%
  left_join(.go)  -> .nodeInfo

.nodeInfo %>% 
  group_by(module) %>%
  summarise(sem(go)) -> .parentTerm

.parentTerm %>% right_join(.nodeInfo) -> .parentTerm

p <- ggplot2::ggplot(.parentTerm, aes(x = V1, y = V2, color = as.factor(cluster))) +
  ggplot2::geom_point(ggplot2::aes(size = pageR), alpha = 0.5) + 
  ggplot2::scale_size_continuous(range = c(0, 10)) +
  ggplot2::scale_x_continuous(name = "") + 
  ggplot2::scale_y_continuous(name = "") + 
  ggplot2::theme_minimal(base_family='GillSans') + 
  ggplot2::theme(legend.position = 'top',
    axis.line.x = ggplot2::element_blank(), 
    axis.line.y = ggplot2::element_blank(),
    strip.background = element_rect(fill = 'grey86', color = 'white')
  ) 
  # scale_color_manual('',values = structure(psME, names = psME))

p <- p + facet_wrap(~module)

data_subset <- .parentTerm %>% ungroup() %>% distinct(parentTerm, .keep_all = T)

p + ggrepel::geom_label_repel(aes(label = parentTerm), 
  data = data_subset,
  # max.overlaps = 5,
  box.padding = grid::unit(1, "lines"), size = 3) +
  labs(x = 'Dimension 1', y = 'Dimension 2') -> p

p + theme(legend.position = 'none')

# nodeInfo %>%
#   filter(go_freq != 1) %>%
#   arrange(desc(go_freq)) %>%
#   head(1) %>%
#   group_by(nodeName) %>%
#   summarise(sem(go)) -> parentTerm

# Grouped go nodes ----

nodeInfo %>%
  filter(go_freq != 1) %>%
  arrange(desc(go_freq)) %>%
  head(1) %>%
  group_by(nodeName) %>%
  summarise(sem(go)) -> data

nodeInfo %>% 
  group_by(module) %>%
    summarise(sem(unlist(go))) -> parentTerm


write_excel_csv(g_df, file = paste0(path, paste(psME, collapse = '_'), '_nodes.csv'))

# Data viz net ----
# 
# g %>% activate('nodes') %>% 
#   mutate(col = ifelse(genus %in% 'Arabidopsis', 'Arabidopsis', domain)) -> g
# 
# g %>%
#   mutate(genus = ifelse(genus %in% ".", '', genus),
#     domain = ifelse(genus %in% ".", '', domain)) -> g

layout = create_layout(g, layout = 'igraph', algorithm = 'kk')

g %>% activate('nodes') %>% distinct(module) %>% pull() -> nodeColor

ggraph(layout) +
  # geom_edge_arc(aes(edge_alpha = weight, color = type), strength = 0.1,) + # edge_width
  geom_node_point(aes(color = col, size = degree)) + 
  geom_node_text(aes(label = genus), repel = TRUE, size = 2) +
  scale_edge_width(range = c(0.3, 1)) +
  theme_graph(base_family = "GillSans") +
  guides(fill=guide_legend(nrow = 2)) +
  theme(legend.position = "top") +
  coord_fixed() +
  scale_color_brewer(type = 'qual')
# scale_color_manual('', values = structure(nodeColor, names = nodeColor) )