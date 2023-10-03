
# current biom file = OTU_level5.biom on Desktop

## 4/20/22
## Microbiome analyses
##


# Start Up----
# import QIIME output into R
#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("biomformat")

#install.packages('biom',repo='http://cran.wustl.edu')

# load packages
library(biomformat)     # used to import QIIME output
library(microbiome)
library(tidyverse)
library(vegan)
library(ggplot2)
library(wesanderson)
library(reshape2)
#library(dplyr)


# import data
  # school desktop
#biom <- read_biom("C:/Users/christine.cornish/Google Drive/PhD/Projects/Ch1_Microbial Microcosms/ch1_data/ch1_16Ssequencing/QIIME2/New/otu_lev5.biom")
  # laptop
biom <- read_biom("/Users/cmc/Library/CloudStorage/GoogleDrive-christine.cornish@ndsu.edu/My Drive/PhD/Projects/Ch1_Microbial Microcosms/OTU_level5.biom")
  # Family level

# extract OTU counts from biom table
otus <- as.matrix(biom_data(biom))

# transpose so that rows are samples and columns are OTUs
otus <- t(otus)

# reorder sample ids -- ascending is default
# otus <- otus[ order(row.names(otus)), ]


# histogram of sample depth
#depths <- rowSums(otus)
#hist(depths, breaks = 30)

# histogram of OTU frequencies
#counts <- colSums(otus > 0)
#hist(counts, breaks = 30)

# remove Unassigned
otus <- otus[ , -652]     # 96 samples, 651 OTUs

# remove OTUs present in <10% of samples
otus <- otus[, colMeans(otus > 0) >= .1]
#depths <- rowSums(otus)
#dim(otus)     # matrix is 96 samples, 431 OTUs

# extract otu file to excel for supplementary table 
otus <- otus %>%
  as.data.frame()

#library(writexl)
#write_xlsx(otus, "/Users/cmc/Library/CloudStorage/GoogleDrive-christine.cornish@ndsu.edu/My Drive/PhD/Projects/Ch1_Microbial Microcosms/otu_family.xlsx")

# replot with removed singletons
#counts <- colSums(otus > 0)
#hist(counts, breaks = 30)

# theme for figures
fig_theme <- function(){
  theme_bw() +
    theme(text = element_text(family = "Times New Roman",
                              color = "black"),
          axis.text = element_text(color = "black",
                                   size = 16),
          axis.title = element_text(size = 16),
          #          axis.line.x = element_line(color="black"),
          #          axis.line.y = element_line(color="black"),
          #          panel.border = element_blank(),
          panel.grid = element_blank(),
          #          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          #          plot.title = element_text(size = 18, vjust = 1, hjust = 0),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16))
  #          legend.position = c(0.95, 0.15),
  #          legend.key = element_blank(),
  #          legend.background = element_rect(color = "black",
  #                                           fill = "transparent",
  #                                           size = 2, linetype = "blank"))
}

# label for facet wrap on figures
day.labs <- as_labeller(
   c("0" = "Pre-treatment", "1" = "2-hours Post-treatment", "14" = "2-weeks Post-treatment"))

trt.labs <- as_labeller(
  c("ctl" = "Control", "gly" = "Glyphosate",
    "frm" = "Commercial Formula", "amp" = "AMPA"))

#
# remove samples with low depth
  # code not necessary as all samples >1000
# sort(depths) [1:10]
# otus <- otus[depths >=1000,]
# dim(otus)
# 

# load metadata
  # school desktop
#map <- read.table("C:/Users/christine.cornish/Google Drive/PhD/Projects/Ch1_Microbial Microcosms/ch1_data/ch1_16Ssequencing/QIIME2/New/metadata.txt",
#                  sep = '\t', head=T, row=1, check=F, comment='') %>%
#  select(-number) %>%
#  mutate(time = as.factor(time),
#         trt = as.factor(trt),
#         lvl = as.factor(lvl))
# laptop
map <- read.table("/Users/cmc/Google Drive (christine.cornish@ndsu.edu)/PhD/Projects/Ch1_Microbial Microcosms/ch1_data/ch1_16Ssequencing/QIIME2/New/metadata.txt",
                  sep = '\t', head=T, row=1, check=F, comment='') %>%
#  select(-number) %>%
  mutate(time = as.factor(time),
         trt = as.factor(trt),
         lvl = as.factor(lvl))

# confirm that metadata and OTU table contain the same samples in the same order
sample.ids <- intersect(rownames(otus), rownames(map)) %>%
  sort()

# subset using sample IDs
otus <- otus[sample.ids,]
map <- map[sample.ids,]
#dim(otus)
#dim(map)


# Species Richness----
rich <- specnumber(otus) %>%
  as.data.frame()

all.equal(rownames(rich), rownames(map))

# merge relative abundance with metadata
rich.meta <- cbind(map, rich)

rich.stats <- rich.meta %>%
  group_by(time, trt, lvl) %>% 
  summarize(mean=mean(.), sd = sd(.))

(box.rich <- rich.meta %>%
    mutate(lvl = recode(lvl, "C" = "0 ppm", "L" = "0.07 ppm", "M" = "0.7 ppm", "H" = "7 ppm"),
           trt = recode(trt, "ctl" = "Control", "gly" = "Glyphosate", "frm" = "Commercial\nFormula", "amp" = "AMPA")) %>%
    mutate(trt = fct_relevel(trt, "Control", "Glyphosate", "Commercial\nFormula",  "AMPA"),
           lvl = fct_relevel(lvl, "0 ppm", "0.07 ppm", "0.7 ppm", "7 ppm")) %>%
    ggplot(aes(x=trt, y=.,
               color = lvl)) +
    geom_boxplot(width = 1.5/length(unique(rich.meta$trt)),
                 position = position_dodge(width = 1, preserve = "total")) +
    geom_jitter(aes(fill = trt), shape = 16, size = 2.5,
                position = position_jitterdodge(dodge.width = 1),
                show.legend = F) +
    geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dashed", size = 0.25) +
    labs(x = "\nTreatment", y = "Observed OTU\n",
         color = "Concentration") +
    scale_color_manual(values = wes_palette(n = 4, name = "Moonrise2")) +
    fig_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(color = "black"),
          strip.text = element_text(size = 14, family = "Times New Roman")) +
    facet_wrap(~time, labeller = day.labs))

ggsave("Fig3.1.jpeg", width = 14, height = 10,
       plot = box.rich)


# Alpha Diversity----
library(tidyverse)
library(ggplot2)
library(vegan)
library(gridExtra)


# calculate shannon H
alpha.div <- otus %>%
  diversity(index = "shannon") %>%
  bind_cols(., map) %>%     # merge diversity calcs with metadata
  rename(shannon = ...1)     # rename diversity calc column

# summary stats of shannon H
alpha.stats <- alpha.div %>%
  group_by(time, trt, lvl) %>% 
  summarize(mean=mean(shannon), sd = sd(shannon))

# tally each trt for n=x in legend label
trt.tally <- alpha.div %>%
  group_by(time, trt) %>%
  summarise(n = n())

# PLOT
# point and errorbar
#(ggpoint.alpha <- alpha.stats %>%
#    mutate(time = recode(time, "0" = "Pre-treatment",
#                         "1" = "2-hours Post-treatment",
#                         "14" = "2-weeks Post-treatment")) %>%
#    mutate(trt = fct_relevel(trt, "ctl", "gly", "amp", "frm")) %>%
#    ggplot(aes(x=time, y=mean)) +
#    geom_errorbar(aes(ymin = mean-sd,
#                      ymax = mean+sd,
#                      color = trt, width = 0.3),
#                  size = 1,
#                  position = position_dodge(width = 1)) +
#    geom_point(aes(color = trt), size = 2,
#               position = position_dodge(width = 1)) +
#    geom_vline(xintercept = c(1.5, 2.5), linetype = "dashed", size = 0.25) +
#    labs(x = "\nTime of Sampling", y = "Average Shannon Diversity\n") +
#    scale_colour_manual(name = "Treatment", labels = c("Control (n=5)", "Glyphosate (n=9)", "AMPA (n=9)", "Formula (n=9)"),
#                        values=wes_palette(n = 4, name = "Moonrise2")) +
#    fig_theme())

#ggsave("alpha_point.png", plot = ggpoint.alpha)

# boxplot
(ggbox.alpha <- alpha.div %>%
    mutate(lvl = recode(lvl, "C" = "0 ppm", "L" = "0.07 ppm", "M" = "0.7 ppm", "H" = "7 ppm"),
           trt = recode(trt, "ctl" = "Control", "gly" = "Glyphosate", "frm" = "Commercial\nFormula", "amp" = "AMPA")) %>%
    mutate(trt = fct_relevel(trt, "Control", "Glyphosate", "Commercial\nFormula",  "AMPA"),
           lvl = fct_relevel(lvl, "0 ppm", "0.07 ppm", "0.7 ppm", "7 ppm")) %>%
    ggplot(aes(x=trt, y=shannon,
               color = lvl)) +
    geom_boxplot(width = 1.5/length(unique(alpha.div$trt)),
                 position = position_dodge(width = 1, preserve = "total")) +
    geom_jitter(aes(fill = trt), shape = 16, size = 2.5,
                position = position_jitterdodge(dodge.width = 1),
                show.legend = F) +
    geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dashed", size = 0.25) +
    labs(x = "\nTreatment", y = "Shannon Diversity\n",
         color = "Concentration") +
    scale_color_manual(values = wes_palette(n = 4, name = "Moonrise2")) +
    fig_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(color = "black"),
          strip.text = element_text(size = 14, family = "Times New Roman")) +
    facet_wrap(~time, labeller = day.labs))

ggsave("Fig2.1.jpeg", width = 14, height = 10,
       plot = ggbox.alpha)


# Relative Abundance----

library(funrar)

# turn into relative abundance
relab.otus <- otus %>%
  as.matrix() %>%
  make_relative()

  # check total abundance in each sample
  apply(relab.otus, 1, sum)
  # check transformed data
  relab.otus[1:5, 1:5]

# check that otus and metadata are in same order
  # output will = TRUE or FALSE
all.equal(rownames(relab.otus), rownames(map))

# merge relative abundance with metadata
relab.meta <- cbind(map, relab.otus) %>%
  rownames_to_column("sample")     # make row headers first column

# reshape data into long format
relab.long <- relab.meta %>% 
  gather(key = "family", value = "relab", 6:436) %>%
  arrange(match(trt, c("ctl", "gly", "amp", "frm"))) #%>%      # sorts columns by trt
  select(-c(time, lvl))

# separate QIIME2 naming into only family name
relab.fam <- relab.long %>%
  separate(col = family, into = c("k", "p", "c", "o", "family"),
           sep = ";") %>%
  select(-c("k", "p", "c", "o"))

# remove values = 0
# no.zeros <- long %>%
#  filter(!relab <= 0) %>%
#  filter(!family == "Unassigned;__;__;__;__") %>%
#  filter(!family == "D_0__Bacteria;__;__;__;__")

# pull top 10 out for legend
leg <- c("Hydrogenophilaceae",
         "Sulfurovaceae",
         "Desulfobacteraceae")
# obtain top ten most abundant taxa
core.taxa <- relab.fam %>%
  group_by(family) %>%
  summarize(sum = sum(relab)) %>%
  top_n(3, sum) %>%
#  pull(family) %>%
  as.data.frame()


# stacked bar plot with all 96 samples
(ggbar.relab <- relab.fam %>%
    mutate(trt = fct_relevel(trt, "ctl", "gly", "frm", "amp")) %>%
    ggplot() +
    geom_bar(aes(x = sample, y = relab, 
                 fill = family),
             position = "stack",
             stat = "identity") +
    labs(x = " ", y = "Relative Abundance\n",
         fill = "Family") +
    scale_y_continuous(expand = c(0,0),
                       limits = c(0,1)) +     # remove gaps in axis
#    scale_fill_discrete(breaks = core.taxa) +     # only show core taxa )
    facet_grid(~trt, scales = "free",
               labeller = trt.labs) +     # scales = free --> removes empty spaces from taxa not shown
    fig_theme() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.text = ggtext::element_markdown(leg),
          strip.background = element_rect(color = "black"),
          strip.text = element_text(size = 14, family = "Times New Roman"),
  legend.position = "bottom"))

ggsave("Fig1.1.jpeg", width = 25, height = 12,
       plot = ggbar.relab)

library(ggh4x)
library(writexl)
write_xlsx(relab.fam, "/Users/cmc/Library/CloudStorage/GoogleDrive-christine.cornish@ndsu.edu/My Drive/PhD/Projects/Ch1_Microbial Microcosms/relab.barplot.xlsx")

# replace "old_value" with "new_value" in column "my_column" using dplyr
relab.fam.1 <- relab.fam %>%
  mutate(family = ifelse(family == "__", "D_4__Unknown Family", family))

relab.fam.2 <- relab.fam.1 %>%
  separate(col = family, into = c("x", "family"),
         sep = "D_4__") %>%
  select(-"x")

# relab plot separated by treatment AND concentration
(ggbar.relab_R2 <- relab.fam.2 %>%
#    filter(!trt == "ctl") %>%
    mutate(lvl = fct_relevel(lvl, "L", "M", "H")) %>%
    mutate(trt = fct_relevel(trt, "gly", "frm", "amp")) %>%
    ggplot() +
    geom_bar(aes(x = sample, y = relab, 
                 fill = family),
             position = "stack",
             stat = "identity") +
    labs(x = " ", y = "Relative Abundance\n",
         fill = "Family") +
    scale_y_continuous(expand = c(0,0),
                       limits = c(0,1)) +     # remove gaps in axis
    #    scale_fill_discrete(breaks = core.taxa) +     # only show core taxa )
    facet_nested_wrap(~ lvl + trt, scales = "free_x", ncol = 3,
               labeller = labeller(trt = trt.labs, 
                                   lvl = c("L" = "0.07 ppm", "M" = "0.7 ppm", "H" = "7 ppm",
                                           "C" = "0 ppm"))) +     # scales = free --> removes empty spaces from taxa not shown
    fig_theme() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
#          legend.text = ggtext::element_markdown(leg),
          strip.background = element_rect(color = "black", fill = "white"),
          panel.border = element_rect(color = "black"),
          strip.text = element_text(size = 14, family = "Times New Roman")))

# extract legend from ggplot
grob <- ggplotGrob(ggbar.relab_R2)
legend_index <- which(grob$layout$name == "guide-box")
legend_grob <- grob$grobs[[legend_index]]

# export legend as image file
png("relab.legend.png", width = 3750, height = 750) # adjust width and height as needed
grid.newpage()
grid.draw(legend_grob)
dev.off() # close the device


ggsave("Fig1_R2.jpeg", width = 45, height = 45,
       plot = ggbar.relab_R2)


# Hierarchical clustering----
# calculate Bray-Curtis distance
otu.dist <- vegdist(relab.otus, method = "bray")

# cluster communities using average-linkage algorithm
otu.clust <- hclust(otu.dist, method = "average")

# plot cluster diagram
plot(otu.clust,
     ylab = "Bray-Curtis dissimilarity")



# Ordination----

# use extracted, transposed otus + load map data

# find overlapping samples in otus and map
common.ids <- intersect(rownames(map), rownames(otus))

# extract all overlapping samples
otus <- otus[common.ids,]
map <- map[common.ids,]

# check that otus and metadata are in same order
# output will = TRUE or FALSE
all.equal(rownames(otus), rownames(map))

# merge relative abundance with metadata
metadata <- cbind(map, otus)

# get Bray-Curtis distances
dist <- vegdist(otus)

# NMDS
set.seed(2)
nmds <- metaMDS(otus,
                k = 2, trymax = 100)     # number of reduced dimensions

# Shepard/stressplot
  # shows scatter around the regression between interpoint distances
    # against their original dissimilarities
stressplot(nmds)

# export information from nmds to plot
# extract NMDS scores (x and y coordinates)
data.scores <- as.data.frame(scores(nmds)$sites)

# add columns from original data to new nmds coordinates
data.scores$time = metadata$time
data.scores$trt = metadata$trt
data.scores$lvl = metadata$lvl

# ordination plot using trt and time
(gg.nmds <- data.scores %>%
    mutate(trt = fct_relevel(trt, "ctl", "gly", "frm", "amp")) %>%
    mutate(time = recode(time, "0" = "Pre-treatment",
                         "1" = "2-hours Post-treatment",
                         "14" = "2-weeks Post-treatment")) %>%
    ggplot(aes(x = NMDS1, y = NMDS2)) +
    geom_point(size = 4, 
               aes(colour = trt, shape = time)) +
    labs(x = "\nNMDS1", y = "NMDS2\n")  +
    guides(colour = guide_legend("Treatment", order = 2),
           shape = guide_legend("Day of\nSampling", order = 1)) +
    scale_colour_manual(labels = c("Control", "Glyphosate", "Commercial Formula", "AMPA"),
                        values=wes_palette(n = 4, name = "Moonrise2")) +
    annotate("text", x = 0.5, y = 0.3, label = "Stress = 0.1793708",
             colour = "black", size = 5, family = "Times New Roman") +
    fig_theme())

# save
ggsave("Fig4.jpeg", width = 11, height = 9,
       plot = gg.nmds)

# ordination plot using trt and lvl, faceted by time
(gg.nmds_R2 <- data.scores %>%
    mutate(trt = fct_relevel(trt, "ctl", "gly", "frm", "amp")) %>%
    mutate(lvl = fct_relevel(lvl, "C", "L", "M", "H")) %>%
    ggplot(aes(x = NMDS1, y = NMDS2)) +
    geom_point(size = 4, 
               aes(colour = trt, shape = lvl)) +
    labs(x = "\nNMDS1", y = "NMDS2\n")  +
    guides(colour = guide_legend("Treatment", order = 2),
           shape = guide_legend("Concentration", order = 1)) +
    scale_colour_manual(labels = c("Control", "Glyphosate", "Commercial\nFormula", "AMPA"),
                        values=wes_palette(n = 4, name = "Moonrise2")) +
    scale_shape_manual(limits = c("C", "L", "M", "H"),
                       labels = c("0 ppm", "0.07 ppm", "0.7 ppm", "7 ppm"),
                       values = c(16, 15, 17, 13)) +
  fig_theme() +
    theme(strip.background = element_rect(color = "black"),
          strip.text = element_text(size = 14, family = "Times New Roman")) +
    facet_grid(~time, labeller = day.labs))

# save
ggsave("SI_Fig2.jpeg", width = 14, height = 10,
       plot = gg.nmds_R2)

# Statistical Analyses----

# run model to check anova assumptions
mod.adiv <- lm(shannon ~ trt*lvl + time, data = alpha.div)
mod.rich <- lm(. ~ trt, data = rich.meta)
mod.rich2 <- lm(. ~ lvl, data = rich.meta)

# normal distribution plot
plot(mod.adiv, which = 2)
plot(mod.rich, which = 2)
plot(mod.rich2, which = 2)

# confirm with shapiro wilk hypothesis test
shapiro.test(alpha.div$shannon)
shapiro.test(rich.meta$.)
  # p is < 0.05, which means these data do NOT come from a normally distributed population

# test out data transformations
log <- log10(alpha.div$shannon)
sqrt <- sqrt(alpha.div$shannon)
cube <- alpha.div$shannon^(1/3)

# non-parametric one-way anova
  # diversity
kruskal.test(shannon ~ trt, data = alpha.div)
kruskal.test(shannon ~ time, data = alpha.div)
kruskal.test(shannon ~ lvl, data = alpha.div)

  # richness
kruskal.test(. ~ trt, data = rich.meta)
kruskal.test(. ~ time, data = rich.meta)
kruskal.test(. ~ lvl, data = rich.meta)

# parametric ANOVA comparison
one.way <- aov(shannon ~ trt, data = alpha.div)
summary(one.way)
two.way <- aov(shannon ~ trt + lvl, data = alpha.div)
summary(two.way)
blocking <- aov(shannon ~ trt + lvl + time, data = alpha.div)
summary(blocking)
interaction <- aov(shannon ~ trt*lvl*time, data = alpha.div)
summary(interaction)

install.packages("AICcmodavg")
library(AICcmodavg)

model.set <- list(one.way, two.way, interaction, blocking)
model.names <- c("one.way", "two.way", "interaction", "blocking")
aictab(model.set, modnames = model.names)
# based on model, blocking is the best fit for data

