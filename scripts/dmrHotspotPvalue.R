#!/gsc/software/linux-x86_64-centos7/R-4.0.2/bin/Rscript --vanilla

#Note: these packages need to be installed.
suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

option_list = list(
  make_option(c("-d", "--density"), type="character", default=NULL, 
              help="density list file", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="sample name", metavar="character"),
  make_option(c("-e", "--event"), type="character", default=NULL, 
              help="event name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="Output directory name", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

d <- read.delim(opt$density,header = T)
sample <- opt$sample
event <- opt$event

ncontrols <- length(d$controls[d$region == "hpv_region"])

de <- d %>%
  group_by(region) %>%
  summarise(sample = c("hpv_sample", rep("other", ncontrols)), 
            dmr.density = c(unique(test), controls))

### ----------------------------------------------------------
### STATS
### ----------------------------------------------------------

fc_mat <- data.frame(region = unique(de$region)) 
fc_mat$log2FC <- NA
fc_mat$zscore <- NA

for (region in unique(de$region)) {
  sub <- de[de$region == region,]
  samp <- sub$dmr.density[sub$sample == "hpv_sample"]
  mean <- mean(sub$dmr.density[sub$sample == "other"])
  fc <- log2(samp / mean)
  z <- (samp-mean(sub$dmr.density))/sd(sub$dmr.density)
  fc_mat[fc_mat$region == region, "log2FC"] <- fc
  fc_mat[fc_mat$region == region, "zscore"] <- z
}

fc <- fc_mat$log2FC[fc_mat$region == "hpv_region"]
num_higher <- nrow(subset(fc_mat, log2FC > fc))
p <- ifelse(num_higher == 0, 0.001, num_higher / 1000)

zs <- fc_mat$zscore[fc_mat$region == "hpv_region"]
num_higher_z <- nrow(subset(fc_mat, zscore > zs))
p2 <- ifelse(num_higher_z == 0, 0.001, num_higher_z / 1000)

pval <- data.frame(sample = sample, event = event, log2FC = fc, pval = p, pval_zscore = p2)

write.table(pval, paste0(opt$out, "/pvalue.txt"), quote = F, col.names = T, sep = "\t", row.names = F)

# add columns and save the density values for plotting later
deSub <- de %>% filter(region == "hpv_region")
#deSub$sample <- sample
deSub$event <- event
write.table(deSub, paste0(opt$out, "/plottingDensityValues.txt"), quote = F, col.names = T, sep = "\t", row.names = F)

### ----------------------------------------------------------
### FIGURE
### ----------------------------------------------------------

plot <- ggplot(de %>% filter(region == "hpv_region"), aes(x = region, y = dmr.density)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(colour = sample), height = 0, width = 0.2, size =3, alpha=0.5) + 
  theme_minimal() +
  labs(x = NULL, y = "DMR Density") +
  annotate("text", x=1.3, y=max(de$dmr.density[de$region == "hpv_region"]), 
           label= paste0("pvalue = ", p), 
           colour = "grey30", size = 5) + 
  scale_colour_manual(values = c("#e63946", "#1d3557")) +
  theme(panel.grid = element_blank(), 
        axis.text.y = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.title = element_text(size=14), 
        axis.ticks.y = element_line(),
        axis.line = element_line(), 
        legend.position = "none")

ggsave(filename = paste0(opt$out, "/regionDensityPlot.pdf"), plot = plot, height = 4, width = 4.5, units = "in")
