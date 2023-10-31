#!/gsc/software/linux-x86_64-centos7/R-4.0.2/bin/Rscript --vanilla
.libPaths("/projects/vporter_prj/R/x86_64-centos7-linux-gnu-library/4.0")

### ----------------------------------------------------------
### OptParse Options
### ----------------------------------------------------------
library(optparse)

option_list = list(
  make_option(c("-d", "--dmr"), type="character", default=NULL, 
              help="DMR file from DSS", metavar="character"),
  make_option(c("-v", "--hpv"), type="character", default=NULL, 
              help="HPV integration regions bed file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="Output directory", metavar="character"),
  make_option(c("-c", "--chrom"), type="character", default=NULL, 
              help="Chromosome sizes", metavar="character"),
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

suppressMessages(library(GenomicRanges))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))

### ----------------------------------------------------------
### Read in Files
### ----------------------------------------------------------

# get the inputs
dmr <- read.delim(opt$dmr, header = T)
out <- opt$out
hpv <- read.delim(opt$hpv, header = F)

# chromosome lengths
cl <- read.delim(opt$chrom, header = F)
sl <- cl$V2
names(sl) <- cl$V1

### ----------------------------------------------------------
### Function
### ----------------------------------------------------------

breakpointHotspotter <- function(breaks, sample){
  
  #Take in all breakpoints and calculate p-values
  gr.list <- GRanges(breaks, seqlengths = sl)
  bw=1000000
  pval=1e-8
  names(gr.list) <- NULL
  
  #gr <- do.call(c, gr.list)
  gr <- GenomicRanges::sort(gr.list)
  
  ## Iterate over chromosomes and calculate p-values
  pranges.list <- GenomicRanges::GRangesList()
  hotspots=data.frame()
  count=1
  
  df = data.frame(chr=c(),midpoint=c(),density=c(),pvalue=c(),null_midpoints=c(),null_density=c())
  
  
  for (chrom in seqlevels(gr)) {
    grc <- gr[seqnames(gr)==chrom]
    
    if (length(grc)>1) {
      nbins <- round(sl[chrom] / 10000)
      midpoints <- (start(grc)+end(grc))/2
      kde <- stats::density(midpoints,bw=bw,kernel='gaussian', from = 0, to = sl[chrom], n = nbins)
      
      # Random distribution of genomic events
      kde.densities <- numeric()
      
      midpoints.r <- round(stats::runif(length(midpoints),1,seqlengths(gr)[chrom]))
      kde.r <- stats::density(midpoints.r,bw=bw,kernel='gaussian', from = 0, to = sl[chrom], n = nbins)
      kde.densities <- c(kde.densities, kde.r$y)
      
      # Use ecdf to calculate p-values 
      p <- 1-stats::ecdf(kde.densities)(kde$y)
      pvalues <- data.frame(chromosome=chrom,start=kde$x,pvalue=p)
      
      # Make GRanges
      pvalues$end <- pvalues$start
      pvalues$chromosome <- factor(pvalues$chromosome, levels=seqlevels(gr))
      pvalues <- as(pvalues,'GRanges')
      seqlevels(pvalues) <- seqlevels(gr)
      
      suppressWarnings(
        seqlengths(pvalues) <- seqlengths(gr)[names(seqlengths(pvalues))]
      )
      
      # Resize from pointsize to bandwidth
      suppressWarnings(
        pvalues <- GenomicRanges::resize(pvalues, width=bw, fix='center')
      )
      pvalues <- trim(pvalues)
      
      ## Find regions where p-value is below specification
      mask <- pvalues$pvalue <= pval
      rle.pvals <- rle(mask)
      rle.pvals$values <- cumsum(rle.pvals$values+1)
      pvalues$group <- inverse.rle(rle.pvals)
      
      if (length(which(mask))>0) {
        pvalues.split <- split(pvalues[mask],pvalues$group[mask])
        pranges <- unlist(endoapply(pvalues.split, function(x) { y <- x[1]; end(y) <- end(x)[length(x)]; y$pvalue <- min(x$pvalue); return(y) }))
        pranges$group <- NULL
        pranges$num.events <- GenomicRanges::countOverlaps(pranges,grc)
        pranges.list[[chrom]] <- pranges
      }
      
      for (el in 1:length(pranges)){
        tmp = as.data.frame(gr[queryHits(findOverlaps(gr,pranges[el],type = "any"))])
        tmp$count=count
        hotspots <- rbind(tmp,hotspots)
        count=count+1
      }
      
      
      mp = kde$x
      ds = kde$y
      pv = pvalues$pvalue
      exp = data.frame(mp=mp,ds=ds)
      exp$type="DMRs"
      exp$p = pv
      
      null_mp = kde.r$x
      null_ds = kde.r$y
      null = data.frame(mp=null_mp,ds=null_ds)
      null$type = "NULL"
      null$p = pv
      
      tmp = rbind(exp,null)
      tmp$chr = chrom
      df <- rbind(tmp,df)
    }
    count=count+1
  }
  pranges <- unlist(pranges.list, use.names=FALSE)
  names(pranges) <- NULL
  
  # make a summarized bed file of the hotspots
  hs <- hotspots %>%
    group_by(count) %>%
    summarise(chr = unique(seqnames), start = min(start), end = max(end))
  
  # write the tables
  write.table(df, paste0(out, "/densityPvalueSummary.txt"), col.names = T,row.names = F, quote = F,sep="\t")
  write.table(hotspots, paste0(out, "/densityDMRHotspots.txt"), col.names = T,row.names = F, quote = F,sep="\t")
  write.table(hs[,2:4], paste0(out, "/densityDMRHotspotsRegions.bed"), col.names = F,row.names = F, quote = F,sep="\t")
  
  return(list(df, hs))
}

### ----------------------------------------------------------
### Run Functions
### ----------------------------------------------------------

bh <- breakpointHotspotter(dmr, out)

### -----------------------------------------------------------------
### PRELIM TEST ON THE REGION
### -----------------------------------------------------------------

b <- dmr
b$midpoint <- round(b$start + ((b$end - b$start)/2))

bedHPV <- b %>%
  filter(chr == hpv$V1[1] & midpoint > hpv$V2[1] & midpoint < hpv$V3[1])

bedChr <- b %>%
  filter(chr == hpv$V1[1])

# chromosome density
d1 <- nrow(bedChr) / cl$V2[cl$V1 == hpv$V1[1]]

# region density
d2 <- nrow(bedHPV) / (hpv$V3 - hpv$V2)

# how much more enriched is the ecDNA vs the chromosome
df <- data.frame(chrDensity = d1, ecDensity = d2)
df$enrichment <- df$ecDensity / df$chrDensity

write.table(df, paste0(out, "/densityDMREnrichment.txt"), col.names = T,row.names = F, quote = F,sep="\t")

### -----------------------------------------------------------------
### PLOTTING THE DMR HOTSPOTS
### -----------------------------------------------------------------

densitySummary <- bh[[1]]
densitySummary <- densitySummary[complete.cases(densitySummary),]

# get the midpoint of the HPV regions
hpv$midpoint <- (hpv$V2 + hpv$V3)/2

# Find the max y value
max_xy <-max(densitySummary$ds)

# change to MB
densitySummary$mp_Mb = densitySummary$mp/1e+06

# Set factor levels
densitySummary$chr <- factor(densitySummary$chr,levels = c(paste0("chr", 1:22),"chrX"))

# line for the HPV integration events
line <- data.frame(xint=(hpv$midpoint/1000000),chr=hpv$V1)
line$chr <- factor(line$chr,levels = c(paste0("chr", 1:22),"chrX"))

# the hotspots bars
hs <- bh[[2]]
hs <- hs[complete.cases(hs),]
hs$chr <- factor(hs$chr,levels = c(paste0("chr", 1:22),"chrX" ))

#### PLOT ALL CHROMOSOMES
densityALL <- ggplot(densitySummary)+
  geom_point(aes(x = mp_Mb, y = ds, color = type),size=1)+
  geom_rect(data = hs, aes(xmin = start/1000000, xmax = end/1000000, ymin = max_xy - (max_xy/15), ymax = max_xy)) +
  facet_wrap(~chr,scales="free")+
  theme_bw() +
  scale_y_continuous(limits = c(0,max_xy)) +
  scale_colour_manual(values = c("#e63946", "#1d3557")) +
  geom_vline(data = line, aes(xintercept = xint), linetype = 2) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=16),
        axis.title =element_text(size=18),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(size=16),
        legend.position = c(0.9, 0.12))+
  labs(y="DENSITY",x="CHROMOSOME POSITION (MB)")+
  guides(color = guide_legend(override.aes = list(size=5)))

ggsave(plot = densityALL, filename = paste0(out,"/densityPlotDMRs.png"), width = 14, height = 8, units = "in")
ggsave(plot = densityALL, filename = paste0(out,"/densityPlotDMRs.pdf"), width = 14, height = 8, units = "in")

### -----------------------------------------------------------------
### SUBSET THE DENSITY AT HPV REGIONS
### -----------------------------------------------------------------

hpv_bins <- list()
for (i in 1:nrow(hpv)) {
  sub <- hpv[i,]
  chrL <- cl$V2[cl$V1 == sub$V1]
  if(sub$midpoint > 3000000 & (chrL - sub$midpoint) > 3000000){
    ds <- densitySummary[densitySummary$chr == sub$V1 & densitySummary$type == "DMRs",]
    ds$dist <- abs(ds$mp - sub$midpoint)
    
    # get the change in density between NULL and DMRs
    #ds <- spread(ds, key = type, value = ds)
    #colnames(ds)[6] <- "null"
    #ds$dsNorm <- ds$DMRs / ds$null
    
    if(nrow(ds) > 0){
      # normalize the density to a position between 0-1 of the max density in the chromosome
      max_ds <- max(ds$ds)
      ds$dsNorm <- ds$ds/max_ds
    
      # find the bin that has HPV integration 
      rownames(ds) <- 1:nrow(ds)
      hpv_bin <- as.numeric(rownames(ds[ds$dist == min(ds$dist),]))
      
      # grab 300 bins (3Mb) on either side of HPV integration
      hpv_df <- ds[(hpv_bin-300):(hpv_bin+300),]
      hpv_df$bin_num <- 1:nrow(hpv_df)
      hpv_df$bin_type <- ifelse(hpv_df$bin_num == 301, "HPV", "flanking")
    }
    else{
      hpv_df <- data.frame(mp=NA, ds=NA, type=NA, p=NA, chr=NA, dist=NA, dsNorm=NA, bin_num=NA, bin_type=NA)
    }
  }
  else{
    hpv_df <- data.frame(mp=NA, ds=NA, type=NA, p=NA, chr=NA, dist=NA, dsNorm=NA, bin_num=NA, bin_type=NA)
  }
  hpv_bins[[i]] <- hpv_df 
}

if(length(hpv_bins) > 0){
  names(hpv_bins) <- hpv$V4
  hpvBins <- bind_rows(hpv_bins, .id = "sites")
  # save the table
  write.table(hpvBins, paste0(out, "/densityHPVBins.txt"), col.names = T,row.names = F, quote = F,sep="\t")  
}


