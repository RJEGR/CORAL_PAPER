

rm(list = ls()) # Limpiar la memoria de la sesion de R
if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE) # 

library(tidyverse)

path <- "~/Downloads/CORAL/"

blastp <- readRDS(paste0(path, 'blastp.rds'))
allData_df <- readRDS(paste0(path,'allData_df.rds'))
load(paste0(path, 'moduleTraitRelationship.RData'))
egg_df <- read_rds(paste0(path, 'trinotate_eggnog.rds'))

# 
# df1 %>% filter(corPvalueStudent < 0.05) %>% pull(module) %>% unique() -> psME
# 
# bwModuleCol %>% filter(module %in% psME) -> pbwModuleCol


ranks <- c()
# Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo

table(blastp$domain)

bwModuleCol %>% 
  left_join(blastp, by = "transcript") %>% 
  filter(domain != '.') -> modules_annot

modules_annot %>% group_by(module) %>% 
  count(domain, sort = T) %>% 
  drop_na(domain) -> alluv_df

library(ggalluvial)

M <- df1 %>% distinct(module) %>% pull

alluv_df %>%
  # arrange(domain) %>%
  ggplot(
  aes(axis1 = module, axis2 = domain, y = n)) +
  geom_alluvium(aes(fill = module), decreasing = T ) +
  geom_stratum(absolute = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
    absolute = FALSE) +
  theme_bw() +
  scale_fill_manual('',values = structure(M, names = M))

genusf <- c("Homo", "Mus", "Danio", "Drosophila", "Arabidopsis")

# modules_annot %>% 
#   filter(domain %in% 'Eukaryota') %>%
#   group_by(module) %>% 
#   count(genus, sort = T) %>% 
#   drop_na(genus) %>% #view()
#   ggplot(aes(x = module, y = n)) +
#   geom_col() +
#   facet_grid(~ genus)

# where Arabidopsis

modules_annot %>% 
  filter(genus %in% genusf) %>%
  group_by(module) %>% 
  count(genus, sort = T) %>% 
  drop_na(genus) # %>% view()


modules_annot %>% 
  # filter(genus %in% c('Arabidopsis', 'Homo')) %>%
  ggplot(aes(y = identity, x = domain)) +
  geom_boxplot()


modules_annot %>% 
  filter(identity > 50) %>%
  filter(domain %in% 'Eukaryota') %>%
  group_by(module) %>% 
  count(genus, sort = T) %>% 
  drop_na(genus) # %>% view()


modules_annot %>% 
  filter(domain %in% 'Eukaryota') %>%
  filter(identity > 50) %>%
  ggplot(aes(y = identity, x = genus)) +
  geom_boxplot() +
  coord_flip()

# by egg

modules_kegg <- bwModuleCol %>% left_join(egg_df) %>% drop_na(Kegg) 

# modules_kegg %>% group_by(module) 

cogs <- read_rds(paste0(path, '/cogs.rds'))

# Process cogs patways

`Information Processing` <- c('J', 'A', 'K', 'L', 'B')
`Cellular Processes` <- c('D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'U', 'O')
`Metabolism` <- c('C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q')

cogsL <- c("Information Processing", "Cellular Processes","Metabolism")

cogs[cogs$code %in% `Information Processing`, 'Pathway'] <- "Information Processing"
cogs[cogs$code %in% `Cellular Processes`, 'Pathway'] <- "Cellular Processes"
cogs[cogs$code %in% `Metabolism`, 'Pathway'] <- "Metabolism"

query.ids <- modules_kegg %>% distinct(transcript) %>% pull(transcript)


mod <- modules_kegg %>% distinct(module) %>% pull()

out <- list()

for(i in mod) {
  j <- i
  
  bwModuleCol %>% filter(module %in% j) %>%
    distinct(transcript) %>%
    pull(transcript) -> which_ids
  
  df <- get_eggnog(egg_df, which_ids) %>% mutate(module = j)
  out[[j]] <- df
  
}

kegg_df <- do.call(rbind, out)


col_palette <- cogs %>% distinct(code, clrs)
col_palette <- structure(col_palette$clrs, names = col_palette$code)


cogs %>% left_join(kegg_df) %>%
  filter(Freq > 5) %>%
  filter(module %in% psME) %>%
  mutate(Pct = Freq/sum(Freq)) %>%
  filter(!grepl('unknown', name)) %>%
  mutate(name = paste(code, name,sep = ', ')) %>%
  mutate(name = fct_reorder(name, Freq)) %>%
  ggplot(aes(y = Freq, x = name)) + 
  geom_col() +
  coord_flip() +
  # scale_y_continuous(labels = scales::percent) +
  labs(x = '' , y = 'Number of transcripts') +
  facet_grid( module ~.) +
  geom_text(aes(label = Freq), size = 2, hjust = -0.05, family = "GillSans") +
  theme_bw(base_family = "GillSans") -> psave

psave + 
  theme(legend.position = 'none',
    strip.background = element_rect(fill = 'grey86', color = 'white'),
    panel.border = element_blank()) -> psave

psave

ggsave(psave, path = path, filename = 'eggnog_bar.png', 
  width = 8, height = 5, dpi = 1000) 
#

# 
# cogs %>% left_join(kegg_df) %>% 
#   arrange(desc(Freq)) %>%
#   mutate(name = factor(name, levels = unique(name))) %>%
#   ggplot(aes(x = name, y = Freq)) +
#   geom_col() +
#   coord_flip() +
#   geom_text(aes(label = Freq), size = 3, hjust = -0.05, family = "GillSans") 

cogs %>% left_join(kegg_df) %>% drop_na(Freq) %>% summarise(sum(Freq))
  # drop_na(Freq) %>%
  filter(Freq > 2) %>% 
  filter(!grepl('unknown', name)) %>% summarise(sum(Freq))
  mutate(name = paste(code, name,sep = ', ')) %>%
  group_by(module) %>%
  mutate(pct = Freq/sum(Freq)) %>%
  mutate(name = fct_reorder(name, Freq)) %>%
  mutate(module = factor(module, levels = hclust$labels[hclust$order])) %>%
  mutate(Pathway = factor(Pathway, levels = cogsL)) %>%
  # filter(module %in% psME) %>%
  # mutate(name = factor(name, levels = unique(name))) %>%
  ggplot(aes(y = name, x = module, fill = pct)) + 
  geom_tile(color = 'white', size = 0.5) +
  geom_text(aes(label = Freq),  vjust = 0.5, hjust = 0.5, size= 4, family = "GillSans") +
  theme_bw(base_family = "GillSans", base_size = 14) +
  ggsci::scale_fill_material(name = "", na.value = "white",
    labels = scales::percent) -> psave

psave <- psave + theme(strip.background = element_rect(fill = 'white', color = 'white'),
  panel.border = element_blank(), axis.text.x = element_text(angle = 20, 
    hjust = 1, vjust = 1, size = 10)) + labs(x = '', y = 'eggNOG')

psave + facet_grid(Pathway ~ ., scales = 'free_y', space = 'free') -> psave

ggsave(psave, path = path, filename = 'eggnog_heatmap.png', 
  width = 9.5, height = 5, dpi = 1000) 

# full eggnog

# kegg_df <- 
kegg_df <- get_eggnog(x = egg_df, ids = unique(egg_df$transcript))


cogs %>% left_join(kegg_df) %>%
  filter(Freq > 2) %>%
  mutate(name = fct_reorder(name, Freq)) %>%
  # filter(!grepl('unknown', name)) %>%
  ggplot(aes(y = Freq, x = name, fill = code)) + 
  geom_col() +
  coord_flip() +
  # scale_y_continuous(labels = scales::percent) +
  labs(x = '' , y = 'Number of transcripts') +
  geom_text(aes(label = Freq), size = 2, hjust = -0.05, family = "GillSans") +
  theme_bw(base_family = "GillSans") -> psave

psave + scale_fill_manual(values = col_palette) +
  theme(legend.position = 'none', panel.border = element_blank()) -> psave

# psave

ggsave(psave, path = path, filename = 'eggnog_global.png', 
  width = 6.5, height = 4.5, dpi = 1000) 

# information storage and processing (J, A, K, L, B), cellular processes and signaling (D, Y, V, T, M, N, Z, U, O), and metabolism (C, G, E, F, H, I, P, Q).

