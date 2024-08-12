#!/gsc/software/linux-x86_64-centos7/R-4.0.2/bin/Rscript --vanilla
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

#cat <- read.delim("/projects/hpv_nanopore_prj/cell_lines/call_integration/output/hela/hpv_size/hpvSizeCategories.txt", header = T)
#out <- "/projects/hpv_nanopore_prj/cell_lines/call_integration/output/hela/hpv_size/hpv_breakpoints.pdf"

# Make help options
option_list = list(
    make_option(c("-c", "--cat"), type="character", default=NULL,
                help="hpvSizeCategories.txt from hpv_integrant_size.R", metavar="character"),
    make_option(c("-o", "--out"), type="character", default = "./",
                help="Output file name", metavar="character")
)

# load in options 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
out <- opt$out
sum <- read.delim(opt$cat, header = T)

# complete integrants
catM <- cat %>% filter(status == "complete")
catM <- separate(catM, col = bp_pair, into = c("bp1", "bp2"), sep = "_", remove = F)
catM <- separate(catM, col = bp1, into = c("bp1_chr", "bp1_pos"), sep = ":", remove = T)
catM <- separate(catM, col = bp2, into = c("bp2_chr", "bp2_pos"), sep = ":", remove = T)

event <- catM
chrs <- unique(c(event$bp1_chr, event$bp2_chr))
sites <- unique(c(paste0(event$bp1_chr,":",event$bp1_pos), paste0(event$bp2_chr,":",event$bp2_pos)))
sitesDF <- data.frame(sites = sites, chr = gsub("\\:.*", "", sites))
    
# dataframes to be filled
intSites <- data.frame(chr = gsub("\\:.*", "", sites), site = as.numeric(gsub(".*:", "", sites)), position = NA)
chrSize <- data.frame(chr = chrs, start = 0, size = NA, intStart = NA, intEnd = NA)
connect <- event[,c("bp1_chr", "bp1_pos", "bp2_chr", "bp2_pos", "category")]
connect$bp1_position <- NA
connect$bp2_position <- NA
    
for (c in chrs) {
    subsite <- sitesDF$sites[sitesDF$chr == c]
        
    # find the minimum position on the chromosome 
    positions <- as.numeric(gsub(paste0(c,":"), "", subsite))
    minsite <- min(positions)
    maxsite <- max(positions)
        
    # add the length to the chrsize
    len <- maxsite - minsite
    chrSize$intStart[chrSize$chr == c] <- minsite
    chrSize$intEnd[chrSize$chr == c] <- maxsite
    chrSize$size[chrSize$chr == c] <- len/1000
        
    # subtract that value from the positions on that chromosome to zero them
    intSites$position[intSites$chr == c] <- (intSites$site[intSites$chr == c] - minsite)/1000
        
    # adjust the positions on the connections
    connect$bp1_position[connect$bp1_chr == c] <- (as.numeric(connect$bp1_pos[connect$bp1_chr == c]) - minsite)/1000
    connect$bp2_position[connect$bp2_chr == c] <- (as.numeric(connect$bp2_pos[connect$bp2_chr == c]) - minsite)/1000
}
    
# count how many connections for each int site
count <- table(c(connect$bp1_pos, connect$bp2_pos))
intSites$nConnect <- count[match(intSites$site, names(count))]
intSites$connect <- ifelse(intSites$nConnect == 1, "one", ifelse(intSites$nConnect == 2, "two", "three_plus"))
    
# seperate heterologous integrants
connectA <- connect[connect$category != "heterologous",]
connectB <- connect[connect$category == "heterologous",]

# title
title <- NULL
for (i in 1:nrow(chrSize)) {
    c <- paste0(chrSize$chr[i], ":", chrSize$intStart[i], "-",chrSize$intEnd[i])
    title[i] <- c
}

    
# make plot
plot <- ggplot() +
    # chromosome bars
    geom_segment(data = chrSize, aes(y = chr, yend = chr, x = start, xend = size), 
                     lineend = "round", color = "lightgrey", size = 3) + 
    geom_point(data = intSites, aes(y = chr, x = position, colour = connect), 
                   size = 3) +
    geom_curve(data = connectA, aes(x = bp1_position, xend = bp2_position, y = bp1_chr, yend = bp2_chr), linetype = 2, size = 0.5, curvature = -0.3, colour = "#fdb714") +
    geom_curve(data = connectB, aes(x = bp1_position, xend = bp2_position, y = bp1_chr, yend = bp2_chr), linetype = 2, size = 0.5, curvature = 0.4, colour = "#219EBC") +
    scale_colour_manual(values = c("black","#BE1E2D", "#ec8992")) +
    theme_bw() +
    labs(x = "position in event (kb)", title = title) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.line = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 10, colour = "black"),
          axis.text.y = element_text(size = 10, colour = "black"),
          axis.text.x = element_text(size = 10, colour = "black"),
          axis.title.x = element_text(size = 12, colour = "black", face = "bold"),
          axis.title.y = element_blank(), legend.position = "none")


ggsave(filename = out, plot = plot, width = 6.5, height = 3.5, units = "in")


