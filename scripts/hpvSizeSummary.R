.libPaths("/projects/vporter_prj/R/x86_64-centos7-linux-gnu-library/4.0")

#Note: these packages need to be installed.
suppressMessages(library(optparse))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(pafr))
suppressMessages(library(tidyverse))
suppressMessages(library(ggsci))
suppressMessages(library(ggpubr))

### -------------------------------------------------------------------------------
### READ IN FILES 
### -------------------------------------------------------------------------------

### CATEGORIES

f <- grep("hpvSizeCategories.txt", c(list.files("/projects/hpv_nanopore_prj/htmcp/call_integration/output", recursive = T, full.names = T),
                                  list.files("/projects/hpv_nanopore_prj/tcga/call_integration/output", recursive = T, full.names = T)), value = T)
n <- gsub("/projects/hpv_nanopore_prj/|tcga|htmcp|/call_integration/output/|/hpv_size/hpvSizeCategories.txt", "", f)

list <- list()
for (i in 1:length(f)) {
    df <- read.delim(f[i], header = T)
    list[[i]] <- df
}
names(list) <- n
cat <- bind_rows(list, .id = "sample")

cat$sample <- gsub("MaM_485", "TCGA-C5-A2LX", cat$sample)
cat$sample <- gsub("MaM_486", "TCGA-C5-A2LY", cat$sample)
int$sample <- gsub("MaM_485", "TCGA-C5-A2LX", int$sample)
int$sample <- gsub("MaM_486", "TCGA-C5-A2LY", int$sample)

### -------------------------------------------------------------------------------
### ADD THE INTEGRATION CATEGORIES
### -------------------------------------------------------------------------------

int <- read.delim("/projects/hpv_nanopore_prj/manuscript2022/tables/integration_types.txt", header = T)
sum <- read.delim("/projects/hpv_nanopore_prj/manuscript2022/tables/intSitesSummaryAll.txt", header = T)

cat$bp1 <- gsub("\\_.*","",cat$bp_pair)
sum$bp <- paste0(sum$chr, ":", sum$pos)

# add the info to the df
cat$event <- sum$event[match(cat$bp1, sum$bp)]
cat$integration_type <- int$type[match(paste0(cat$sample, cat$event), paste0(int$sample, int$event))]
cat$integration_type[is.na(cat$integration_type)] <- "unmatched_integration"
cat$HPV.type <- sum$HPVchr[match(paste0(cat$sample, cat$event), paste0(sum$Patient, sum$event))]
    
### -------------------------------------------------------------------------------
### GET STATS ON THE INTEGRANTS
### -------------------------------------------------------------------------------

# Longest integrant
head(cat %>% arrange(desc(max_nHPV)), 1)
head(cat %>% filter(status == "complete") %>% arrange(desc(max_nHPV)), 1)

# complete integrants
catC <- cat %>% filter(status == "complete")

# percent heterologous
table(catC$category)
table(catC$category)/sum(table(catC$category))

# heterologous 
catCH <- catC %>% filter(category == "heterologous")
table(catCH$integration_type)
table(catCH$HPV.type)
table(catCH$integration_type)/sum(table(catCH$integration_type))

# length per HPV type
catC %>%
    group_by(HPV.type) %>%
    summarise(mean = mean(max_nHPV))

### -------------------------------------------------------------------------------
### FIGURES
### -------------------------------------------------------------------------------

ann_colors <- list(HPV.clade = c(A7="#B55F8F", A9="#253083", Other="#A9A9A9"),
                   HPV.type = c(HPV16="#3953A4", HPV18="#9768ad", HPV45="#CC138C", HPV82="#369797",
                                HPV52="#0B8DCD", HPV31="#d2ecf9", HPV73="#737474", HPV68="#f3b2d4", HPV97="black", HPV58="#8bafba", HPV59="#ae3030"))

# hpv type
ggplot(catC %>% filter(HPV.type %in% c("HPV16", "HPV18")), aes(x = HPV.type, colour = HPV.type, y = max_nHPV))+
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(height = 0, width = 0.2, size =3, alpha=0.5) +
    theme_minimal() +
    xlab("HPV type") + 
    ylab("max # of HPV genomes in integrant") +
    scale_colour_manual(values = c(ann_colors[["HPV.type"]])) +
    theme(panel.grid = element_blank(), 
          axis.text = element_text(size = 12,colour = "black"),
          axis.title = element_text(size=12, face = "bold", colour = "black"), 
          axis.ticks.y = element_line(),
          axis.line = element_line(), 
          legend.position = "none") +
    stat_compare_means(method = "wilcox")

catC$category_simple <- ifelse(catC$category %in% c("partial", "full"), "single", "heterologous")

# integration type
ggplot(catC %>% filter(!integration_type %in% c("complex_ecDNA_integration","unmatched_integration") & category_simple == "heterologous"), aes(x = integration_type))+
    geom_bar() + 
    theme_minimal() +
    xlab("integration type") + 
    ylab("# of heterologous integrants") +
    #scale_colour_manual(values = c(ann_colors[["HPV.type"]])) +
    theme(panel.grid = element_blank(), 
          axis.text.y = element_text(size = 12,colour = "black"),
          axis.text.x = element_text(size = 12,colour = "black", angle = 60, hjust = 1, vjust = 1),
          axis.title = element_text(size=12, face = "bold", colour = "black"), 
          axis.ticks.y = element_line(),
          axis.line = element_line())

catC$all <- "all"
# integration type
ggplot(catC %>% filter(HPV.type %in% c("HPV16", "HPV18")), aes(x = all, fill = category_simple))+
    geom_bar(position="fill", colour = "black") + 
    theme_minimal() +
    xlab("integration type") + 
    ylab("% of integrants") +
    labs(fill = NULL) +
    scale_fill_manual(values = c("#219ebc","#ffb703")) +
    theme(panel.grid = element_blank(), 
          axis.text.y = element_text(size = 12,colour = "black"),
          axis.text.x = element_text(size = 12,colour = "black", angle = 60, hjust = 1, vjust = 1),
          legend.text = element_text(size = 12,colour = "black"),
          axis.title = element_text(size=12, face = "bold", colour = "black"), 
          axis.ticks.y = element_line(),
          axis.line = element_line())


# hpv type
ggplot(catC %>% filter(!integration_type %in% c("complex_ecDNA_integration","unmatched_integration")), aes(x = integration_type, y = max_nHPV))+
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(height = 0, width = 0.2, size =2, alpha=0.5, colour = "dark grey") +
    theme_minimal() +
    xlab("integration type") + 
    ylab("max # of HPV genomes in integrant") +
    #scale_colour_manual(values = c(ann_colors[["HPV.type"]])) +
    theme(panel.grid = element_blank(), 
          axis.text.y = element_text(size = 12,colour = "black"),
          axis.text.x = element_text(size = 12,colour = "black", angle = 60, hjust = 1, vjust = 1),
          axis.title = element_text(size=12, face = "bold", colour = "black"), 
          axis.ticks.y = element_line(),
          axis.line = element_line(), 
          legend.position = "none") 

# ngroups 

ggplot(catC, aes(x = as.factor(ngroups), fill = category_simple))+
    geom_bar(colour = "black") + 
    theme_minimal() +
    xlab("# of integrant structures") + 
    ylab("# of breakpoint pairs") +
    labs(fill = NULL) +
    scale_fill_manual(values = c("#219ebc","#ffb703")) +
    theme(panel.grid = element_blank(), 
          axis.text.y = element_text(size = 12,colour = "black"),
          axis.text.x = element_text(size = 12,colour = "black"),
          legend.text = element_text(size = 12,colour = "black"),
          axis.title = element_text(size=12, face = "bold", colour = "black"), 
          axis.ticks.y = element_line(),
          axis.line = element_line())

# incomplete 

catI <- cat[cat$status == "incomplete",]

ggplot(catI, aes(x = max_length))+
    geom_histogram(colour = "black", size = 0.5) + 
    theme_minimal() +
    xlab("max length of incomplete HPV integrant (bp)") + 
    ylab("count") +
    labs(fill = NULL) +
    theme(panel.grid.minor = element_blank(), 
          axis.text.y = element_text(size = 12,colour = "black"),
          axis.text.x = element_text(size = 12,colour = "black"),
          legend.text = element_text(size = 12,colour = "black"),
          axis.title = element_text(size=12, face = "bold", colour = "black"), 
          axis.ticks.y = element_line(),
          axis.ticks.x = element_line(),
          axis.line = element_line())

    