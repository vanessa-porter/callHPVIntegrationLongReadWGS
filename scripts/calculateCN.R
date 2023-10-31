#!/gsc/software/linux-x86_64-centos7/R-4.0.2/bin/Rscript --vanilla

#Note: these packages need to be installed.
suppressMessages(require(data.table))
suppressMessages(require(optparse))
suppressMessages(require(dplyr))


# Make help options
option_list = list(
  make_option(c("-p", "--ploidy"), type="character", default=NULL,
              help="Ploidy file from PloiDetect", metavar="character"),
  make_option(c("-c", "--cna"), type="character", default=NULL,
              help="CNA segments from PloiDetect", metavar="character"),
  make_option(c("-d", "--depth"), type="character",
              help="Region depth bed file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default = "regionCN.txt",
              help="Output file name", metavar="character")
)

# load in options 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
out <- opt$out

# ploidy file
ploidy_file = fread(opt$ploidy)
ploidy = ploidy_file$ploidy[1]

# segment file
segments = fread(opt$cna)

# get depth at ploidy segments
depth_ploidy_segments = weighted.mean(segments[CN == ploidy]$segment_depth, segments[CN == ploidy]$end - segments[CN == ploidy]$pos)
tc = ploidy_file$tp[1]

# calculate depth of contaminating normal 
normal_dep = depth_ploidy_segments * (2*(1-tc))/(2*(1-tc) + ploidy * tc)

# the difference between ploidy depth and normal depth
diff = (depth_ploidy_segments - normal_dep) / ploidy
regress_deps = seq(from = normal_dep, by = diff, length.out = 10)

# calculate CN at specified regions
## copy number = (depth - normal_dep) / diff

experimental_depths = fread(opt$depth)
names(experimental_depths) = c("chr", "pos", "end", "depth")
experimental_depths[,cn:=(depth - normal_dep)/diff]
experimental_depths$cn <- ifelse(experimental_depths$cn < 0, 0, experimental_depths$cn)
experimental_depths$cn_rounded <- round(experimental_depths$cn/0.5*0.5)

# categorize the regions
outside <- experimental_depths[grep("chr", experimental_depths$chr),]
outside <- experimental_depths$pos[c(1,nrow(outside))]
experimental_depths$region <- NA
experimental_depths$region[experimental_depths$pos %in% outside] <- "outside"
experimental_depths$region[!experimental_depths$pos %in% outside] <- "inside"

# detect which parts are amplified
cn_outside <- mean(experimental_depths$cn_rounded[experimental_depths$region == "outside"])
depth_outside <- mean(experimental_depths$depth[experimental_depths$region == "outside"])

experimental_depths <- experimental_depths %>%
  mutate(
    cn_change = case_when(
      region == "inside" ~ cn_rounded - cn_outside,
      TRUE ~ 0)) %>%
  mutate(
    cn_change_category = case_when(
      region == "inside" & cn_rounded == cn_outside ~ "no_change",
      region == "inside" & cn < cn_outside ~ "loss",
      region == "inside" & cn > cn_outside ~ "gain",
      TRUE ~ "outside"
  )) %>%
  mutate(
    perc_change = case_when(
      region == "inside" ~ depth / depth_outside,
      TRUE ~ 0)) 

# save file
write.table(experimental_depths, file = out, quote = F, col.names = T, row.names = F, sep = "\t")

