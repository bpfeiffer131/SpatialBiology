library(umap)
library(tidyverse)
library(RColorBrewer)
library(gtsummary)
library(patchwork)


# Theme -------------------------------------------------------------------


## Set graphics theme
theme_set(theme_bw() +
            theme(panel.grid = element_blank()))

theme_update(
  axis.text.x = element_text(color = "black", face = "bold", size = 10),
  axis.text.y = element_text(color = "black", face = "bold", size = 10),
  axis.title = element_text(face = "bold"),
  plot.title = element_text(face = "bold"),
  strip.text = element_text(face = "bold"))


# Data load ---------------------------------------------------------------


space.dataset <- readxl::read_xlsx("/Users/brentpfeiffer/Lab/20240101_SpatialBiology/Subset of BM_GVHD_Measurements2.xlsx")

space.dataset <- as.data.frame(unclass(space.dataset),                     # Convert all columns to factor
                             stringsAsFactors = TRUE)

# Concatenate Parent factors
space.dataset <- space.dataset %>% 
  mutate(Group = fct_collapse(Parent, 
                              BM = c("BM_1", "BM_2"),
                              GVHD = c("GVHD_1", "GVHD_2")))

space.data <- space.dataset %>% 
  select(11:45, 47)

space.label <- space.dataset %>% 
  select(1:10, 51)

space.label <- as.data.frame(unclass(space.label),                     # Convert all columns to factor
                       stringsAsFactors = TRUE)

space.label2 <- space.dataset %>% 
  select(1:10, 17, 47:50)


# General Stats -----------------------------------------------------------


# General stats of the dataset

hist(space.dataset$Nucleus..ki67.mean)

space.dataset %>% 
  filter(Cell.Categories %in% c("CD4+ Effector", "CD8+ Effector", 
                                "CD38+ Enterocyte", "CD44+ Glycoprotein", 
                                "Dentritic Cell", "Hematopietic Cell", 
                                "Marcrophage", "Treg", "Unclassified Cell")
           ) %>% 
  mutate(Cell.Categories = fct_infreq(Cell.Categories)) %>% 
  count(Cell.Categories)

space.dataset %>% 
  filter(Cell.Categories %in% c("CD4+ Effector", "CD8+ Effector", 
                                "CD38+ Enterocyte", "CD44+ Glycoprotein", 
                                "Dentritic Cell", "Hematopietic Cell", 
                                "Marcrophage", "Treg", "Unclassified Cell")
  ) %>% 
  mutate(Cell.Categories = fct_infreq(Cell.Categories)) %>% 
  ggplot(aes(y = Cell.Categories, fill = Parent)) +
  geom_bar()
  

# Standard Graphics -------------------------------------------------------


# Treg Ki67 gradient by Nucleus circularity vs groups

tregki <- space.dataset %>% 
  # filter(Nucleus..ki67.mean > 10) %>% 
  filter(Cell.Categories == "Treg") %>% 
  ggplot(aes(Group, Nucleus..Circularity)) +
  geom_violin() +
  geom_boxplot(width = 0.25,
               outlier.shape = NA) +
  # geom_jitter(aes(color = Nucleus..ki67.mean),
  #             width = 0.2,
  #             alpha = 0.7,
  #             size = 2,
  #             stroke = "black") +
  geom_point(aes(color = Nucleus..ki67.mean),
    position = position_jitter(width = .2, seed = 0),
    size = 5, alpha = 1
  ) +
  geom_point(aes(color = Nucleus..ki67.mean),
    position = position_jitter(width = .2, seed = 0),
    size = 5, stroke = .7, shape = 1, color = "black"
  ) +
  scale_color_gradient(low = "grey90", high = "grey10") +
  labs(title = "Treg Nucleus Circularity & Ki67 Expression",
       color = "Mean Ki67",
       y = "Nucleus Circularity") +
  theme(legend.position = c(.9, .2))


ggsave("TregKi67.tiff", units="in", width=10, height=6, dpi=300, compression = 'lzw')

tregki


tregki67 <- space.dataset %>% 
  # filter(Nucleus..ki67.mean > 10) %>% 
  filter(Cell.Categories == "Treg") %>% 
  ggplot(aes(Group, Nucleus..ki67.mean)) +
  geom_hline(yintercept = 10) +
  #geom_violin() +
  #geom_boxplot(width = 0.25,
               #outlier.shape = NA) +
  # geom_jitter(aes(color = Nucleus..ki67.mean),
  #             width = 0.2,
  #             alpha = 0.7,
  #             size = 2,
  #             stroke = "black") +
  geom_point(aes(color = Nucleus..ki67.mean),
             position = position_jitter(width = .3, seed = 0),
             size = 5, alpha = 1
  ) +
  geom_point(aes(color = Nucleus..ki67.mean),
             position = position_jitter(width = .3, seed = 0),
             size = 5, stroke = .7, shape = 1, color = "black"
  ) +
  scale_color_gradient(low = "grey90", high = "grey10") +
  labs(title = "Treg Mean Ki67 Expression vs Group",
       color = "Mean Ki67",
       y = "Mean Nuclear Ki67 Expression") +
  theme(legend.position = c(.9, .7),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave("TregKi67-1.tiff", units="in", width=10, height=6, dpi=300, compression = 'lzw')

tregki67


# CD4 Effector Ki67 gradient by Nucleus circularity vs groups
cd4ki <- space.dataset %>% 
  filter(Nucleus..ki67.mean > 20) %>% 
  filter(Cell.Categories == "CD4+ Effector") %>% 
  ggplot(aes(Group, Nucleus..Circularity)) +
  geom_violin() +
  geom_boxplot(width = 0.25,
               outlier.shape = NA) +
  geom_jitter(aes(color = Nucleus..ki67.mean),
              width = 0.4,
              alpha = 0.7,
              size = 2) +
  scale_color_gradient(low = "grey90", high = "grey10") +
  labs(title = "CD4+ Effector Nucleus Circularity & Ki67 Expression > 20",
       color = "Mean Ki67")

cd8ki <- space.dataset %>% 
  filter(Nucleus..ki67.mean > 20) %>% 
  filter(Cell.Categories == "CD8+ Effector") %>% 
  ggplot(aes(Group, Nucleus..Circularity)) +
  geom_violin() +
  geom_boxplot(width = 0.25,
               outlier.shape = NA) +
  geom_jitter(aes(color = Nucleus..ki67.mean),
              width = 0.4,
              alpha = 0.7,
              size = 2) +
  scale_color_gradient(low = "grey90", high = "grey10") +
  labs(title = "CD8+ Effector Nucleus Circularity & Ki67 Expression > 20",
       color = "Mean Ki67")

macroki <- space.dataset %>% 
  filter(Nucleus..ki67.mean > 10) %>% 
  filter(Cell.Categories == "Macrophage") %>% 
  ggplot(aes(Group, Nucleus..Circularity)) +
  geom_violin() +
  geom_boxplot(width = 0.25,
               outlier.shape = NA) +
  geom_jitter(aes(color = Nucleus..ki67.mean),
              width = 0.4,
              alpha = 0.7,
              size = 2) +
  scale_color_gradient(low = "grey90", high = "grey10") +
  labs(title = "Macrophage Nucleus Circularity & Ki67 Expression > 20",
       color = "Mean Ki67")

(tregki | cd4ki) / (cd8ki | macroki)


space.dataset %>% 
  # filter(Nucleus..ki67.mean > 10) %>% 
  filter(Cell.Categories == "Treg") %>% 
  ggplot(aes(Nucleus..ki67.mean, fill = Group)) +
  geom_histogram(bins = 25) +
  scale_fill_brewer(palette = "Dark2")
  

# Centroid Plot of the whole slide
space.dataset %>% 
  filter(Cell.Categories == "Treg" |
           Cell.Categories == "CD4+ Effector" |
           Cell.Categories == "CD8+ Effector" |
           Cell.Categories == "Macrophage" |
           Cell.Categories == "Dentritic Cell" |
           Cell.Categories == "CD38+ Enterocyte" |
           Cell.Categories == "CD44+ Glycoprotein" |
           Cell.Categories == "Unclassified Cell") %>%
  mutate(Cell.Categories = fct_rev(Cell.Categories)) %>% 
  ggplot(aes(Centroid.X.µm, Centroid.Y.µm)) +
  geom_point(aes(color = Cell.Categories),
             alpha = 0.7, shape = 19) +
  facet_wrap(~ Parent, scales = "free") +
  scale_colour_brewer(palette = "Spectral")


space.dataset %>% 
  filter(Cell.Categories == "Treg" |
           Cell.Categories == "CD4+ Effector" |
           Cell.Categories == "CD8+ Effector" |
           Cell.Categories == "Macrophage" |
           Cell.Categories == "Dentritic Cell") %>%
  mutate(Cell.Categories = fct_rev(Cell.Categories)) %>% 
  ggplot(aes(Centroid.X.µm, Centroid.Y.µm)) +
  geom_point(aes(color = Cell.Categories),
             alpha = 0.7, shape = 19) +
  facet_wrap(~ Parent, scales = "free") +
  scale_colour_brewer(palette = "Spectral")


# UMAP Analysis -----------------------------------------------------------



# Default UMAP - KNN 15 and min-dist 0.1
space.umap <- umap(space.data)

umap_df <- space.umap$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2")

umapdf_combined = dplyr::bind_cols(space.label, umap_df)


umapgraph1 <- umapdf_combined %>% 
  filter(Cell.Categories == "CD8+ Effector" |
           Cell.Categories == "CD4+ Effector" |
           Cell.Categories == "Macrophage" |
           Cell.Categories == "Dentritic Cell" |
           Cell.Categories == "Hematopietic Cell" |
           Cell.Categories == "CD38+ Enterocyte" |
           Cell.Categories == "Treg" |
           Cell.Categories == "CD44+ Glycoprotein") %>%
  mutate(Cell.Categories = fct_relevel(Cell.Categories, 
                                       c("Treg", 
                                        "CD4+ Effector",
                                        "CD8+ Effector",
                                        "Dentritic Cell",
                                        "Macrophage",
                                        "Hematopietic Cell",
                                        "CD44+ Glycoprotein",
                                        "CD38+ Enterocyte"))) %>% 
  ggplot(aes(UMAP1, UMAP2)) +
  geom_point(aes(color = Cell.Categories), alpha = 0.7, size = 2) +
  facet_wrap(vars(Group)) +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "All Cell Categories",
       color = "Cell Categories")

ggsave("UMAP_Plot4.png", units="in", width=10, height=6, dpi=300)

 
umapgraph2 <- umapdf_combined %>% 
  filter(Cell.Categories == "CD8+ Effector" |
           Cell.Categories == "CD4+ Effector" |
           Cell.Categories == "Macrophage" |
           Cell.Categories == "Dentritic Cell" |
           Cell.Categories == "Hematopietic Cell" |
           Cell.Categories == "Treg") %>%
  mutate(Cell.Categories = fct_relevel(Cell.Categories, 
                                       c("Treg", 
                                         "CD4+ Effector",
                                         "CD8+ Effector",
                                         "Dentritic Cell",
                                         "Macrophage",
                                         "Hematopietic Cell"
                                         ))) %>% 
  ggplot(aes(UMAP1, UMAP2)) +
  geom_point(aes(color = Cell.Categories), alpha = 0.7, size = 2) +
  facet_wrap(vars(Group)) +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "Immune Cell Categories",
       color = "Cell Categories") +
  theme(legend.key.width = unit(0.75, 'cm'),
        legend.key.height = unit(1.5, 'cm'),
        legend.title = element_text(size = 14, 'bold'),
        legend.text = element_text(size = 12)
        )

ggsave("UMAP_Plot3.png", units="in", width=10, height=6, dpi=300)

umapgraph1

umapgraph2

umapgraph1 / umapgraph2


# Ki67 information with UMAP

umapdf_combined3 = dplyr::bind_cols(space.label2, umap_df)

Ki67count <- umapdf_combined3 %>% 
  filter(Cell.Categories == "CD8+ Effector" |
           Cell.Categories == "CD4+ Effector" |
           Cell.Categories == "Macrophage" |
           Cell.Categories == "Dentritic Cell" |
           Cell.Categories == "Hematopietic Cell" |
           Cell.Categories == "Treg") %>% 
  mutate(Ki67.Experession = cut(Nucleus..ki67.mean,
                                breaks = c(0, 20, 255),
                                include.lowest = T,
                                right = F))

Ki67count %>% 
  mutate(Group = fct_collapse(Parent, 
                              BM = c("BM_1", "BM_2"),
                              GVHD = c("GVHD_1", "GVHD_2"))) %>% 
  count(Group, Cell.Categories, Ki67.Experession) %>%
  group_by(Group, Cell.Categories) %>% 
  mutate(percent = n / sum(n)) %>% 
  as.matrix() %>% 
  as.table()

Ki67count %>% 
  mutate(Group = fct_collapse(Parent, 
                              BM = c("BM_1", "BM_2"),
                              GVHD = c("GVHD_1", "GVHD_2"))) %>%
  count(Group, Cell.Categories, Ki67.Experession) %>%
  group_by(Group, Cell.Categories) %>% 
  mutate(percent = n / sum(n)) %>% 
  ggplot(aes(Ki67.Experession, Cell.Categories)) +
  geom_tile(aes(fill = percent)) +
  facet_wrap(~Group) +
  scale_fill_stepsn(n.breaks = 10, colours = viridis::viridis(9)) +
  scale_x_discrete(labels = c("[0,20)" = "No Expression", 
                   "[20,255]" = "Mod - High Expression")) +
  labs(title = "Ki67 Expression Percentage by Cell Type",
       fill = "Percent (%)",
       y = "",
       x = "")

Ki67count %>% 
  mutate(Group = fct_collapse(Parent, 
                              BM = c("BM_1", "BM_2"),
                              GVHD = c("GVHD_1", "GVHD_2"))) %>%
  count(Group, Cell.Categories, Ki67.Experession) %>%
  group_by(Group, Cell.Categories) %>% 
  mutate(percent = n / sum(n)) %>% 
  ggplot(aes(Ki67.Experession, Cell.Categories)) +
  geom_tile(aes(fill = percent)) +
  geom_text(aes(label = round(percent, 2))) +
  facet_wrap(~Group) +
  scale_fill_gradient(low = "white", high = "red") +
  #scale_fill_stepsn(n.breaks = 10, colours = viridis::viridis(9)) +
  scale_x_discrete(labels = c("[0,20)" = "No Expression", 
                              "[20,255]" = "Mod - High Expression")) +
  labs(title = "Ki67 Expression Percentage by Cell Type",
       fill = "Percent (%)",
       y = "",
       x = "")

ggsave("Ki67Exp.png", units="in", width=10, height=6, dpi=300)


Ki67count %>% 
  count(Group, Cell.Categories, Ki67.Experession) %>%
  group_by(Group, Cell.Categories) %>% 
  mutate(percent = n / sum(n)) %>% 
  ggplot(aes(Ki67.Experession, Cell.Categories)) +
  geom_tile(aes(fill = percent)) +
  geom_text(aes(label = round(percent, 2))) +
  facet_wrap(~Group) +
  scale_fill_gradient(low = "white", high = "red") +
  #scale_fill_stepsn(n.breaks = 10, colours = viridis::viridis(9)) +
  scale_x_discrete(labels = c("[0,20)" = "No Expression", 
                              "[20,255]" = "Mod - High Expression")) +
  labs(title = "Ki67 Expression Percentage by Cell Type",
       fill = "Percent (%)",
       y = "",
       x = "")


Ki67count %>% 
  select(Cell.Categories, Group, Ki67.Experession) %>% 
  tbl_strata(strata = Group,
             ~.x %>% 
  tbl_summary(by = Ki67.Experession ))

umapdf_combined3 %>% 
  filter(Cell.Categories == "CD8+ Effector" |
           Cell.Categories == "CD4+ Effector" |
           Cell.Categories == "Macrophage" |
           Cell.Categories == "Dentritic Cell" |
           #Cell.Categories == "Hematopietic Cell" |
           Cell.Categories == "Treg" & 
           Nucleus..ki67.mean > 30) %>%
  mutate(Cell.Categories = fct_relevel(Cell.Categories, 
                                       c("Treg", 
                                         "CD4+ Effector",
                                         "CD8+ Effector",
                                         "Dentritic Cell",
                                         "Macrophage",
                                         "Hematopietic Cell"))) %>% 
  ggplot(aes(UMAP1, UMAP2)) +
  geom_point(aes(color = Cell.Categories), alpha = 0.7, size = 2) +
  facet_wrap(vars(Parent)) +
  scale_colour_brewer(palette = "Set1")



# Custom UMAP KNN 15 and min-dist 0.5

custom.config = umap.defaults # Set of configurations
custom.config$min_dist = 0.1 # change min_dist
custom.config$knn

space.umap1 <- umap(space.data, config=custom.config)

umap1_df <- space.umap1$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2")

umapdf_combined1 = dplyr::bind_cols(space.label, umap1_df)


umapdf_combined1 %>% 
  filter(Cell.Categories == "Treg" |
           Cell.Categories == "CD4+ Effector" |
           Cell.Categories == "CD8+ Effector" |
           Cell.Categories == "Macrophage" |
           Cell.Categories == "Dentritic Cell") %>% 
  mutate(Cell.Categories = fct_rev(Cell.Categories)) %>% 
  ggplot(aes(UMAP1, UMAP2)) +
  geom_point(aes(color = Cell.Categories), alpha = 0.7, size = 2) +
  facet_wrap(vars(Parent)) +
  scale_colour_brewer(palette = "Spectral")

umapdf_combined1

# Custom UMAP KNN 100 and min-dist 0.1

custom.config = umap.defaults # Set of configurations
custom.config$n_neighbors = 100 # change KNN

space.umap2 <- umap(space.data, config=custom.config)

umap2_df <- space.umap2$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2")

umapdf_combined2 = dplyr::bind_cols(space.label, umap2_df)


umapdf_combined2 %>% 
  filter(Cell.Categories == "Treg" |
           Cell.Categories == "CD4+ Effector" |
           Cell.Categories == "CD8+ Effector" |
           Cell.Categories == "Macrophage" |
           Cell.Categories == "Dentritic Cell") %>%
  mutate(Cell.Categories = fct_rev(Cell.Categories)) %>% 
  ggplot(aes(UMAP1, UMAP2)) +
  geom_point(aes(color = Cell.Categories), alpha = 0.7, size = 2) +
  facet_wrap(vars(Parent)) +
  scale_colour_brewer(palette = "Spectral")


umapdf_combined2 %>% 
  mutate(Tx.Group = fct_recode(Parent, BM = "BM_1", BM = "BM_2", GHVD = "GVHD_1", 
                            GHVD = "GVHD_2")) %>% 
  filter(Cell.Categories == "Treg" |
           Cell.Categories == "CD4+ Effector" |
           Cell.Categories == "CD8+ Effector" |
           Cell.Categories == "Macrophage" |
           Cell.Categories == "Dentritic Cell") %>%
  ggplot(aes(Cell.Categories)) +
  geom_bar(aes(fill = Tx.Group)) +
  coord_flip() +
  scale_fill_brewer(palette = "Paired")
  #facet_wrap(vars(Tx.Group))

umapdf_combined %>%
  select(Parent, Cell.Categories) %>%
  filter(Cell.Categories == "Treg" |
           Cell.Categories == "CD4+ Effector" |
           Cell.Categories == "CD8+ Effector" |
           Cell.Categories == "Macrophage" |
           Cell.Categories == "Dentritic Cell") %>%
  tbl_summary(by = Parent)
