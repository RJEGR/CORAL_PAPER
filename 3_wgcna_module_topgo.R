


rm(list = ls());

if(!is.null(dev.list())) dev.off()

options(verbose = F)

library(tidyverse)
library(topGO)

path <- "~/Downloads/CORAL/"

datExpr <- readRDS(paste0(path, 'datExpr.rds'))

# make a function to loop per module the enrichment
load(paste0(path, 'moduleTraitRelationship.RData'))

source(paste0(path, 'functions.R'))

allData_df <- readRDS(paste0(path,'allData_df.rds'))

# set significant modules and genes

df1 %>% filter(corPvalueStudent < 0.05) %>% pull(module) %>% unique() -> psME

bwModuleCol %>% filter(module %in% psME) -> pbwModuleCol


# TRINO=~/Documents/Tools/Trinotate-Trinotate-v3.2.1/util
# export PATH=$TRINO:$PATH 
# extract_GO_assignments_from_Trinotate_xls.pl  --Trinotate_xls Pocillopora_annotation_report.xls  -T > Trinotate_report.xls.gene_ontology

go_file <- paste0(path, 'Trinotate_report.xls.gene_ontology')

MAP <- topGO::readMappings(go_file)

out <- list()

# modules <- psME
M <- df1 %>% distinct(module) %>% pull

for(i in M) {
  # # df needs to include follow cols: "ids"            "sampleB"        "GO.ID" (list)      and    "padj" 
  j <- i
  
  df <- allData_df %>% filter(module %in% j)
  
  cat("\nNumber of genes in the module",j, "\n")
  cat(df %>% count() %>% pull())
  cat("\n")
  
  
  allRes <- GOenrichment(df, Nodes = 20)
  
  allRes <- allRes %>% mutate(module = j) %>% as_tibble()
  out[[j]] <- allRes
  
}

# la categoria biological process fue la usada para describir en el art. de Oscar.
# Functions of animal-type gene products were analyzed using animal models as background (human, mouse, zebrafish, and fly), whereas for the functions of plant-type gene products (belonging to the photosynthetic endosymbiont) the Arabidopsis thaliana model was used as background.

do.call(rbind, out) -> topGO

topGO %>% mutate_at(vars(-one_of(c('GO.ID', 'Term', 'module'))), as.double) -> topGO

table(topGO$module) # Top 

# 
# topGO %>%
#   group_by(module) %>%
#   mutate(Size = Annotated / max(Annotated)) %>%
#   mutate(Term = fct_reorder(Term, Size, .desc = F)) %>%
#   # mutate(Term = factor(Term)) %>%
#   ggplot(aes(y = Term, x = Size)) +
#   geom_col() + facet_grid( module ~ ., scales = 'free')

write.csv(topGO, file = paste0(path, 'topGO.csv'))
# Reduce terms 

# Term reduction by

library(GOSemSim)

library(rrvgo)

# reduce_out <- list()

# for(i in psME) {
#   
#   j <- i
#   
#   df <- topGO %>% filter(module == j)
#  
#   allRes <- allRes %>% mutate(module = j) %>% as_tibble()
#   reduce_out[[j]] <- allRes
#   
# }

str(goid <- sort(unique(topGO$GO.ID))) # using only significant genes


scores <- -log(as.numeric(topGO$classicKS))

names(scores) <- topGO$GO.ID

# if arabitopsis org.At.tair.db

SimMatrix <- calculateSimMatrix(goid, orgdb = 'org.Hs.eg.db', ont="BP", method = 'Wang')

reducedTerms <- reduceSimMatrix(SimMatrix, threshold = 0.9, orgdb = 'org.Hs.eg.db')

# and plot
#
reducedTerms %>% 
  dplyr::select(go, cluster, parent, parentTerm) %>% 
  left_join(topGO, by = c('go' = 'GO.ID')) -> reducedTerms


write.csv(reducedTerms, file = paste0(path, 'topGO.csv'))

reducedTerms %>% 
  group_by(module, parentTerm) %>%
  summarise(Annotated = sum(Annotated)) %>%
  group_by(module) %>%
  mutate(Size = Annotated / max(Annotated)) %>%
  mutate(parentTerm = fct_reorder(parentTerm, Size, .desc = T)) %>%
  mutate(module = factor(module, levels = rev(hclust$labels[hclust$order]))) %>%
  ggplot(aes(y = parentTerm, x = module, size = Size)) +
  geom_point() + 
  theme_bw(base_family = "GillSans") -> psave

  
filename <- paste0('reduced_term_', paste(M, collapse = '_'), '.png')

ggsave(psave, filename = filename, path = path, width = 9, height = 5) 

# cmd

y <- cmdscale(as.matrix(as.dist(1 - SimMatrix)), eig = TRUE, k = 2)

data <- cbind(as.data.frame(y$points), 
  reducedTerms[match(rownames(y$points), reducedTerms$go),]) %>%
  # group_by(module, parentTerm) %>%
  group_by(module) %>%
  mutate(Size = Annotated / max(Annotated))


p <- ggplot2::ggplot(data, ggplot2::aes(x = V1, y = V2)) + # color = parentTerm
  ggplot2::geom_point(ggplot2::aes(size = Size, color = module), alpha = 0.5) + 
  ggplot2::scale_size_continuous(range = c(0, 10)) +
  ggplot2::scale_x_continuous(name = "") + 
  ggplot2::scale_y_continuous(name = "") + 
  ggplot2::theme_minimal(base_family='GillSans') + 
  ggplot2::theme(legend.position = 'top',
    axis.line.x = ggplot2::element_blank(), axis.line.y = ggplot2::element_blank()
    ) +
  scale_color_manual('',values = structure(M, names = M))

p + facet_wrap(~parentTerm)
# scale_size(name = '', guide = FALSE, range=c(0, 5))

# p + facet_wrap(~module) -> p

data_subset <- data %>% ungroup() %>% distinct(parentTerm, .keep_all = T)

p + ggrepel::geom_label_repel(aes(label = parentTerm), 
  data = data_subset,
  # max.overlaps = 5,
  box.padding = grid::unit(1, "lines"), size = 3) +
  labs(x = 'Dimension 1', y = 'Dimension 2') -> p


ggsave(p, filename = 'reducedTerms.png', path = path, 
  width =10, height = 9, dpi = 300) 

data %>% group_by(module) %>% count(parentTerm)

