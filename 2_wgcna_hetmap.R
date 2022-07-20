
library(WGCNA)
library(flashClust)
library(tidyverse)


rm(list = ls());

if(!is.null(dev.list())) dev.off()

path <- "~/Downloads/CORAL/"

datExpr <- readRDS(paste0(path, 'datExpr.rds'))

load(paste0(path, 'count_annot_Boquita_vs_Carrizales_up_down_genes.Rdata'))

bwnet <- readRDS(paste0(path, "bwnet.rds"))

bwmodules = labels2colors(bwnet$colors)

names(bwmodules) <- names(bwnet$colors)

table(bwmodules)

# reads <- colSums(datExpr)

reads <- rowSums(count0[rownames(count0) %in% names(bwmodules),])
Total <- sum(reads)

data.frame(reads, bwmodules) %>% 
  as_tibble() %>% 
    group_by(bwmodules) %>% 
  summarise(n = n(), reads = sum(reads)) %>%
  dplyr::rename('module' = 'bwmodules') -> stats

# stats %>% view()

Brown <- grepl('brown', rownames(datExpr))
Green <- grepl('green', rownames(datExpr))
Boquita <- grepl('Boquita', rownames(datExpr))
Carrizales <- grepl('Carrizales', rownames(datExpr))

# datTraits <- data.frame(Brown, Green, Boquita, Carrizales) # esta parece la mejor combinacion de grupos!!!

# if use only site 
datTraits <- data.frame(Boquita, Carrizales)

datTraits <- 1*datTraits

# datTraits <- diag(nrow = 4,ncol = 4)
# colnames(datTraits) <- rownames(datExpr)
rownames(datTraits) <- rownames(datExpr)

# Recalculate MEs with color labels

MEs0 = moduleEigengenes(datExpr, bwmodules)$eigengenes

MEs = orderMEs(MEs0)

# rownames(MEs) <- str_replace_all(rownames(MEs), '^ME', '')
names(MEs) <- str_replace_all(names(MEs), '^ME', '')

moduleTraitCor = cor(MEs, datTraits, use= "p")

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(datTraits))


moduleTraitCor %>% as_tibble(rownames = 'module') %>% 
  pivot_longer(-module, values_to = 'moduleTraitCor') -> df1

moduleTraitPvalue %>% as_tibble(rownames = 'module') %>% 
  pivot_longer(-module, values_to = 'corPvalueStudent') %>%
  right_join(df1) -> df1

hclust <- hclust(dist(moduleTraitCor), "complete")

up_df %>% distinct(transcript) %>% pull() -> updegs
down_df %>% distinct(transcript) %>% pull() -> dwndegs

bwmodules %>% 
  as_tibble(., rownames = 'transcript') %>%
  mutate(degs = ifelse(transcript %in% updegs, 'up', 
    ifelse(transcript %in% dwndegs, 'down', ''))) %>%
  dplyr::rename('module' = 'value') -> bwModuleCol

bwModuleCol %>%  group_by(module) %>% count(sort = T) %>% left_join(stats)

bwModuleCol %>% 
  group_by(module, degs) %>% count(sort = T) -> bwModuleDF

bwModuleDF %>% mutate(module = factor(module, levels = hclust$labels[hclust$order])) -> bwModuleDF

bwModuleDF %>% group_by(module) %>% mutate(pct = n / sum(n)) -> bwModuleDF



df1 %>%
  mutate(star = ifelse(corPvalueStudent <.001, "***", 
    ifelse(corPvalueStudent <.01, "**",
      ifelse(corPvalueStudent <.05, "*", "")))) -> df1

df1 %>%
  mutate(name = factor(name, levels = c('Boquita', 'Carrizales', 'Green', 'Brown'))) %>%
  mutate(facet = ifelse(name %in% c('Boquita', 'Carrizales'), 'Site', 'Morphotype')) %>%
  mutate(moduleTraitCor = round(moduleTraitCor, 2)) %>%
  # mutate(star = ifelse(star != '', paste0(moduleTraitCor, '(', star,')'), moduleTraitCor)) %>%
  mutate(star = ifelse(star != '', paste0(moduleTraitCor, '(', star,')'), '')) %>%
  ggplot(aes(y = module, x = name, fill = moduleTraitCor)) +
  # geom_tile(color = 'black', size = 0.5, width = 0.7) + 
  geom_raster() +
  geom_text(aes(label = star),  vjust = 0.5, hjust = 0.5, size= 5, family =  "GillSans") +
  ggsci::scale_fill_gsea(name = "", reverse = T, na.value = "white") +
  # scale_fill_viridis_c(name = "Membership", na.value = "white") +
  ggh4x::scale_y_dendrogram(hclust = hclust) +
  labs(x = '', y = 'Module') +
  guides(fill = guide_colorbar(barwidth = unit(3, "in"),
    alignd = 0.5,
    ticks.colour = "black", ticks.linewidth = 0.5,
    frame.colour = "black", frame.linewidth = 0.5,
    label.theme = element_text(size = 12))) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'white', color = 'white'),
    axis.line.y = element_line(color = 'white'),
    axis.text.y = element_text(hjust = 1.2),
    axis.ticks.length = unit(5, "pt")) +
    facet_wrap(~ name, scales = 'free_x') -> p1 

p1 <- p1 + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank())

p1 <- p1 + theme(panel.spacing.x = unit(0, "mm"))

ggsave(p1, filename = 'ModuleTraitRelationship_1.png', 
  path = path, width = 5, height = 5, dpi = 1000)


# negative lFC == Carrizales
# Positive lFC == Boquita


bwModuleDF %>% 
  # mutate(degs = ifelse(degs == 'ns', '', degs)) %>%
  ggplot(aes(x = module, y = n, fill = degs)) + 
  labs(y = 'Number of transcripts') +
  geom_col() + coord_flip() +
  # geom_col(color = 'black', size = 0.25) + coord_flip() +
  scale_fill_manual(name = '', values = c("grey30", "#EE4141", "#2428D0")) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top",
    axis.title.y = element_blank(), axis.text.y= element_blank(),
    axis.ticks.y=element_blank(), axis.line.y = element_blank()) -> p2

library(patchwork)

p1 + p2 + plot_layout(widths = c(1, 0.5)) +
  labs(caption = '* corPvalueStudent < 0.05 ') -> psave

ggsave(psave, filename = 'ModuleTraitRelationship_2.png', path = path, width = 6.7, height = 5,dpi = 1000)

# set significant modules

df1 %>% filter(corPvalueStudent < 0.05) %>% pull(module) %>% unique() -> psME

bwModuleCol %>% filter(module %in% psME) -> pbwModuleCol

# get relevant genes

head(bwmodules_sig <- bwmodules[bwmodules %in% psME])

save(moduleTraitCor, moduleTraitPvalue, MEs, psME, datTraits, df1, bwModuleDF, bwModuleCol,hclust,
  file = paste0(path, 'moduleTraitRelationship.RData'))


# Intramodular connectivity ----
# We would find modules containing genes w/ high positive / negative correlations in spite of a variable (for example morphotype awa site) previosly correlated with plotEigengeneNetworks

# Calculate the correlations between modules

geneModuleMembership <- as.data.frame(WGCNA::cor(datExpr, MEs, use = "p"))

# pca of modules

PCA <- prcomp(t(geneModuleMembership), scale. = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
# sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dtvis <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])


dtvis %>%
  mutate(id = rownames(.)) %>%
  ggplot(., aes(PC1, PC2)) +
  # geom_point(size = 5, alpha = 0.9) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  geom_text(aes(label = id), alpha = 0.9) +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) 
  # theme_bw(base_family = "GillSans", base_size = 16) +
  # theme(plot.title = element_text(hjust = 0.5), legend.position = 'top')


# What are the p-values for each correlation?
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)))

# What's the correlation for the trait?
geneTraitSignificance <- as.data.frame(cor(datExpr,datTraits, use = "p"))

# What are the p-values for each correlation?
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples = nrow(datExpr)))
names(geneTraitSignificance) <- paste("GS.", names(datTraits), sep = "")
names(GSPvalue) <- paste("p.GS.", names(datTraits), sep = "")


GSpval <- GSPvalue %>%
  tibble::rownames_to_column(var = "ids")

gMM_df <- geneModuleMembership %>%
  tibble::rownames_to_column(var = "ids") %>%
  gather(key = "module", value = "moduleMemberCorr", -ids) 

gMM_df %>% 
  mutate(module = str_replace_all(module, '^ME', '')) %>%
  # mutate(module = factor(module, levels = yLabels)) %>% 
  as_tibble() -> gMM_df

# Prepare gene significance df
GS_bacprod_df <- geneTraitSignificance %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "ids")


# Put everything together 
allData_df <- gMM_df %>%
  left_join(GS_bacprod_df, by = "ids") %>%
  left_join(GSpval, by = "ids") %>%
  as_tibble()

# Filter by significance genes 

# allData_df %>% filter(ids %in% names(bwmodules_sig)) -> allData_df

bwmodules %>% as_tibble(rownames = 'ids') %>% rename('module' = 'value') %>%
  left_join(allData_df) -> allData_df

saveRDS(allData_df, file = paste0(path,'allData_df.rds'))

# Sanity check
allData_df %>% group_by(module) %>% count() 

# Longer
allData_df %>% select_at(vars(starts_with('GS'))) %>% names() -> colName

allData_df %>%
  select_at(vars(starts_with('GS'), 'module')) %>%
  # mutate(n = 1:nrow(.)) %>%
  pivot_longer(cols = all_of(colName), names_to = 'group', values_to = 'GS') %>%
  mutate(group = str_replace_all(group, '^GS.', '')) %>% # pairs = paste0(module,'-', group)
  mutate(GS = abs(GS)) -> gs_df

allData_df %>% select_at(vars(starts_with('p.GS'))) %>% names() -> colName

allData_df %>%
  select_at(vars(starts_with('p.GS'), 'moduleMemberCorr')) %>%
  # mutate(n = 1:nrow(.)) %>%
  pivot_longer(cols = colName, names_to = 'group', values_to = 'p.GS') %>%
  mutate(group = str_replace_all(group, '^p.GS.', '')) %>% 
  select(p.GS, moduleMemberCorr) %>%
  cbind(gs_df, .) -> df_viz

plotscatterSig(df_viz, M = psME)


allData_df %>%
  ggplot(aes(moduleMemberCorr)) +
  geom_histogram() +
  facet_wrap(~ module)
