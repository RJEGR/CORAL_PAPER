
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

# g %>% activate("edges") %>% filter(!is.na(weight)) 


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
    # walktrap = cluster_walktrap(.)$membership,
    # fast_greedy = cluster_fast_greedy(.)
    ) -> g

hist(V(g)$louvain)
# hist(V(g)$walktrap)
table(V(g)$domain)

g %>% activate("nodes") %>% as_tibble() -> nodeInfo


write_rds(g, file = paste0(path, paste(psME, collapse = '_'), '_network.rds'))

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


sem <- function(x, top = T) {
  # if(n > 1) {
  require(rrvgo)
  
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
    
    if(top) {
      data %>%
        arrange(desc(size)) %>%
        summarise(parent_freq = length(unique(parent)),
          parentTerm = head(.,1)$parentTerm,
          size = head(.,1)$size,
          across(parent, .fns = llist)) -> data 
    } else {
      data <- data
    }
      

    
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

  return(data)
  
}


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
  # head(1) %>%
  group_by(nodeName) %>%
  summarise(sem(go)) -> parentTerm

# file <- paste0(path, paste(psME, collapse = '_'), '_parentTerm.rds')

# write_rds(parentTerm, .parentTerm, file = file)


.parentTerm %>% ungroup() %>% 
  select(nodeName, parentTerm, size, parent)  %>% 
  group_by(nodeName) %>%
  summarise(parent_freq = length(unique(parent)),
    parentTerm = parentTerm,
    size = size,
    across(parent, .fns = llist)) %>%
  select(names(parentTerm)) -> .parentTerm

rbind(parentTerm, .parentTerm) %>%
  # distinct(nodeName, parentTerm) %>%
  left_join(nodeInfo) -> df

drop_go <- function(x) {paste(x, sep = ';', collapse = ';')}

df %>%
  mutate(parent = drop_go(parent), 
         go = drop_go(go)) -> df

write_excel_csv(df, file = paste0(path, paste(psME, collapse = '_'), '_reduceSimGO_nodes.csv'))

# Continue w/  -----

df <- read_csv(paste0(path, paste(psME, collapse = '_'), '_reduceSimGO_nodes.csv'))

df %>% count(module)

# 1 blue        132
# 2 turquoise   170

df %>% distinct(parentTerm) # todavia se puede reducir terminos parent !!!

x <- df$parent
x <- sort(strsplit(unlist(x), ";")[[1]])

hsGO <- read_rds(paste0(path, '/hsGO_BP.rds'))

SimMatrix <- calculateSimMatrix(x, orgdb = 'org.Hs.eg.db', ont="BP", 
  semdata = hsGO, method = 'Wang')

reducedTerms <- reduceSimMatrix(SimMatrix, threshold = 0.9, orgdb = 'org.Hs.eg.db')

y <- cmdscale(as.matrix(as.dist(1 - SimMatrix)), eig = TRUE, k = 2)

data <- cbind(as.data.frame(y$points), 
  reducedTerms[match(rownames(y$points), reducedTerms$go),])

data 

# 
# df %>% 
#   mutate(parent = list(parent)) %>%
#   group_by(module) %>%
#   summarise(sem(parent, top = FALSE)) -> .df


df %>% ggplot(aes(size, degree)) + geom_point() + 
  facet_wrap(~ module, scales = 'free_y')

df %>%
  filter(size > 0) %>%
  distinct(parentTerm, .keep_all = T) %>%
  group_by(module) %>%
  slice_max(n = 20, order_by = size) %>%
  # arrange(size) %>%
  # mutate(parentTerm = factor(parentTerm, levels = parentTerm)) %>%
  mutate(parentTerm = fct_reorder(parentTerm, pageR, .desc = F)) %>%
  ggplot(aes(y = pageR, x = parentTerm)) +
  geom_col() +
  coord_flip() +
  ggsci::scale_fill_material() +
  labs(x = 'GO parentTerm (top 20)', y = 'Centrality degree') +
  facet_grid(module ~ ., scales = 'free') +
  theme_classic(base_family = "GillSans") +
  theme(strip.background = element_rect(fill = 'grey86', color = 'white'),
    panel.border = element_blank(), legend.position = 'top') -> p

p

filename <- paste0(paste(psME, collapse = '_'), '_reduceSimGO_nodes.png')

ggsave(p, path = path, filename = filename, 
  width = 5, height = 5, dpi = 1000) 

# Single batch ----

g %>% activate("nodes") %>% as_tibble() -> nodeInfo

nodeInfo %>%
  group_by(module) %>%
  drop_na(go) %>%
    summarise(sem(unlist(go))) -> data

data %>% count(module)
# 
# p <- ggplot2::ggplot(data, aes(x = V1, y = V2, 
#   color = as.factor(cluster))) +
#   ggplot2::geom_point(ggplot2::aes(size = size), alpha = 0.5) + 
#   ggplot2::scale_size_continuous(range = c(0, 10)) +
#   # ggplot2::scale_x_continuous(name = "") + 
#   # ggplot2::scale_y_continuous(name = "") + 
#   ggplot2::theme_bw(base_family='GillSans', base_size = 14) + 
#   ggplot2::theme(legend.position = 'top',
#     axis.line.x = ggplot2::element_blank(), 
#     axis.line.y = ggplot2::element_blank(),
#     strip.background = element_rect(fill = 'grey86', color = 'white'),
#     panel.border = element_blank()
#   ) 
# 
# p <- p + facet_wrap(~module)
# 
# data_subset <- data %>% distinct(parentTerm, .keep_all = T)
# 
# p + ggrepel::geom_text_repel(aes(label = parentTerm), 
#   data = data_subset,
#   family = 'GillSans',
#   max.overlaps = 10,
#   box.padding = grid::unit(1, "lines"), size = 5) +
#   labs(x = 'Dimension 1', y = 'Dimension 2') -> p
# 
# p + theme(legend.position = 'none') + ggsci::scale_color_jco() -> psave
# 
# filename <- paste0(paste(psME, collapse = '_'), '_reduceSimGO.png')
# 
# ggsave(psave, path = path, filename = filename, 
#   width = 12, height = 10, dpi = 1000) # 
# 
# data_subset %>%
#   # filter(size > 0) %>%
#   group_by(module) %>%
#   arrange(desc(size)) %>%
#   mutate(parentTerm = fct_reorder(parentTerm, size, .desc = F)) %>%
#   ggplot(aes(y = size, x = parentTerm)) +
#   geom_col() +
#   coord_flip() +
#   labs(x = 'GO parentTerm (top 20)') +
#   facet_grid(module ~ ., scales = 'free') +
#   theme_classic(base_family = "GillSans") +
#   theme(strip.background = element_rect(fill = 'grey86', color = 'white'),
#     panel.border = element_blank(), legend.position = 'top') -> p
# 
# filename <- paste0(paste(psME, collapse = '_'), '_reduceSimGO_nodes.png')
# 
# ggsave(p, path = path, filename = filename, 
#   width = 5, height = 5, dpi = 1000) 

# Data viz net ----

# g %>% activate("nodes") %>% as_tibble() -> nodeInfo
# 
# nodeInfo %>%
#   group_by(module) %>%
#   drop_na(go) %>%
#   summarise(sem(unlist(go))) -> df

df <- read_csv(paste0(path, paste(psME, collapse = '_'), '_reduceSimGO_nodes.csv'))
g <- read_rds(paste0(path, paste(psME, collapse = '_'), '_network.rds'))

g %>% activate('nodes') %>% left_join(df %>% select(nodeName, parentTerm)) -> g

g %>% activate('edges') %>% filter(!is.na(weight)) -> g

g %>% activate('nodes') %>% as_tibble() %>% count(module)

g %>%
  activate('nodes') %>%
  mutate(hub = centrality_hub(),
    degree = centrality_degree(),
    eigen = centrality_eigen(),
    pageR = centrality_pagerank()) -> g

hist(V(g)$pageR)

degree_f <- function(x) {quantile(x, probs = 0.95)}


g %>% activate('nodes') %>% 
  group_by(module) %>%
  mutate(nodeName = ifelse(hub > degree_f(hub), nodeName, NA)) %>%
  as_tibble() %>% drop_na(nodeName) %>% count(module) 

g %>% activate('nodes') %>% 
  group_by(module) %>%
  mutate(nodeName = ifelse(hub > degree_f(hub), nodeName, NA)) %>%
  ungroup() %>%
  filter(!is.na(nodeName)) -> g

layout = create_layout(g, layout = 'igraph', algorithm = 'kk')
# layout = create_layout(g, layout = 'igraph', algorithm = 'star')

# g %>% activate('nodes') %>% distinct(module) %>% pull() -> nodeColor


ggraph(layout) +
  facet_nodes(~ module) +
  scale_color_manual('', values = structure(psME, names = psME) ) +
  geom_edge_arc(aes(edge_alpha = weight), strength = 0.1) + # edge_width
  geom_node_point(aes(color = module, size = pageR)) +
  # geom_node_text(aes(label = parentTerm), repel = TRUE, size = 2) +
  ggrepel::geom_text_repel(data = layout, aes(x, y, label = parentTerm))
  scale_edge_width(range = c(0.3, 1)) +
  theme_graph(base_family = "GillSans") +
  guides(fill=guide_legend(nrow = 2)) +
  theme(legend.position = "top") +
  coord_fixed() +
  scale_color_brewer(type = 'qual')
# scale_color_manual('', values = structure(nodeColor, names = nodeColor) )
