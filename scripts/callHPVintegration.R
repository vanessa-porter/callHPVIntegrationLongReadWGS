#!/gsc/software/linux-x86_64-centos7/R-4.0.2/bin/Rscript --vanilla

#Note: these packages need to be installed.
suppressMessages(library(optparse))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(pafr))
suppressMessages(library(bedtoolsr))
options(bedtools.path = "/path/to/bedtools-2.27.1/bin")

# Make help options
option_list = list(
  make_option(c("-v", "--vcf"), type="character", default=NULL,
              help="Sniffles VCF file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="hpv_out",
              help="Output directory name", metavar="character")
)

# load in required 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# read in necessary files
vcf <- read.delim(opt$vcf, comment.char = '#', header = F, sep = "\t", stringsAsFactors = F)

# make directory
outdir <- opt$outdir
dir.create(outdir) # makes in current directory 

# Extract info on the SV calls
info <- strsplit(vcf$V8, ";")

# extract the reads and make a list with each SV call
all_reads <- list()
for (i in 1:length(info)){
  read <- gsub(".*=", "", info[[i]][grep("RNAMES", info[[i]])])
  all_reads[[i]] <- unlist(strsplit(read, ","))
}

# Add the SV type and chr2 pos2 to the vcf dataframe
vcf$SV.type <- gsub(".*=", "", lapply(info, function(x){x[grep("SVTYPE",x)]}))
vcf$chr2 <- gsub(".*=", "", lapply(info, function(x){x[grep("CHR2",x)]}))
vcf$end <- gsub(".*=", "", lapply(info, function(x){x[grep("END",x)]}))
vcf$VAF <- gsub(".*=", "", lapply(info, function(x){x[grep("AF=",x)]}))

# Only include the conventional chromosomes
vcf <- vcf[vcf$V1 %in% c(paste0("chr",as.character(1:22)), "chrX"),]

# add an ID number for the SV calls
vcf$SV.id <- paste0(vcf$SV.type, rownames(vcf))

# subset the translocations
tra <- vcf[vcf$SV.type == "TRA",]

# Make contingincies if there are 0, 1, and >1 HPV integration sites in the sample
if (nrow(tra[grep("HPV",tra$chr2),]) == 0)  {
  
  print("This sample has zero HPV integration sites with 5+ reads")
  
  events <- data.frame(chr=NA, pos=NA, HPVchr=NA, HPVpos=NA, 
                       HPV.site="none", event="none", VAF=NA, read.depth=0)
  
  # save summary file
  write.table(events, file = paste0(outdir,"/summary.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
  
} else  if (nrow(tra[grep("HPV",tra$chr2),]) == 1) {
  
  print("This sample has 1 HPV integration site with 5+ reads")
  
  hpv_sites <- tra[grep("HPV",tra$chr2),]
  
  # get the number of reads in the site
  hpv_info <- unlist(strsplit(hpv_sites$V8, ";"))
  read <- gsub(".*=", "", hpv_info[10])
  hpv_reads <- unlist(strsplit(read, ","))
  
  # Make a summary file
  events <- data.frame(chr=hpv_sites$V1, pos=hpv_sites$V2, HPVchr=hpv_sites$chr2, HPVpos=hpv_sites$end, 
                       HPV.site="hpv.site1", event="event1", VAF=hpv_sites$VAF, read.depth=length(hpv_reads))
  
  # save summary file
  write.table(events, file = paste0(outdir,"/summary.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
  
  # HPV-integration site reads
  read <- gsub(".*=", "", strsplit(hpv_sites$V8, ";")[[1]][10])
  hpv_reads <- unlist(strsplit(read, ","))
  write(hpv_reads, file = paste0(outdir,"/hpv.site1.txt"), sep = "\t")
  
  # integration site bed file
  events_bed <- events[,c(1,2,2,5)]
  events_bed$pos.1 <- events_bed$pos.1 + 1
  write.table(events_bed, file = paste0(outdir,"/hpv_integration_sites.bed"), sep = "\t", quote = F, col.names = F, row.names = F)
  
  # Find associated SVs
  event_SV <- filter(vcf, grepl(paste(hpv_reads, collapse="|"), V8))
  write.table(event_SV[,1:10], file = paste0(outdir,"/hpv.site1.vcf"), sep = "\t", quote = F, col.names = F, row.names = F)
  
} else {
  
  # take out the HPV integration sites
  hpv_sites <- tra[grep("HPV",tra$chr2),]
  rownames(hpv_sites) <- 1:nrow(hpv_sites)
  hpv_sites$SV.id <- paste0("hpv.site",rownames(hpv_sites))
  
  # HPV-integration site reads
  hpv_info <- strsplit(hpv_sites$V8, ";")
  
  hpv_reads <- list()
  for (i in 1:length(hpv_info)){
    read <- gsub(".*=", "", hpv_info[[i]][10])
    hpv_reads[[i]] <- unlist(strsplit(read, ","))
  }
  names(hpv_reads) <- hpv_sites$SV.id
  
  ####
  #### REMOVE DUPLICATE SITES (within 3 bp)
  ####
  hpv_sites <- arrange(hpv_sites, V1, V2)
  
  collapse_sites <- function(hpv_sites, hpv_reads){
    diff <- diff(hpv_sites$V2, lag = 1)
    repeats <- which(diff < 4 & diff >= 0)
    
    # get the sites needed to be merged
    repeats <- as.list(repeats)
    sites <- lapply(repeats, function(x){
      x <- c(x,x+1)
      x <- hpv_sites$SV.id[x]
    })
    
    # merge the reads from them
    new_reads <- lapply(sites, function(x){
      x <- Reduce(c,hpv_reads[x])
    })
    
    # select the site with the most reads to be the main break site
    kept_site <- lapply(sites, function(x){
      a <- length(unlist(hpv_reads[x[1]]))
      b <- length(unlist(hpv_reads[x[2]]))
      ifelse(a > b, x[1], x[2])
    })
    names(new_reads) <- Reduce(c,kept_site)
    
    # remove the old sites
    remove <- Reduce(c,sites)[!Reduce(c,sites) %in% Reduce(c,kept_site)]
    keep <- hpv_sites$SV.id[!hpv_sites$SV.id %in% remove]
    hpv_reads_filt <- hpv_reads[keep]
    
    # replace the contents 
    for (site in Reduce(c,kept_site)) {
      hpv_reads_filt[[site]] <- new_reads[[site]]
    }
    hpv_sites <- hpv_sites[hpv_sites$SV.id %in% keep,]
    
    return(list(hpv_reads_filt, hpv_sites))
  }
  
  filt <- collapse_sites(hpv_sites = hpv_sites, hpv_reads = hpv_reads)
  filt2 <- collapse_sites(hpv_sites = filt[[2]], hpv_reads = filt[[1]])
  filt3 <- collapse_sites(hpv_sites = filt2[[2]], hpv_reads = filt2[[1]])
  filt4 <- collapse_sites(hpv_sites = filt3[[2]], hpv_reads = filt3[[1]])
  
  hpv_reads_filt <- filt4[[1]]
  hpv_sites <- filt4[[2]]
  
  print(paste0("This sample has ", length(hpv_reads_filt), " HPV integration sites with 5+ reads"))
  
  ####
  #### MAKE A READ MATRIX FOR HPV INTEGRATION SITES
  ####
  
  # Make a matrix that shows which events each read spans
  mat <- matrix(nrow = length(unique(as.character(unlist(hpv_reads_filt)))), ncol = nrow(hpv_sites), dimnames = list(unique(as.character(unlist(hpv_reads_filt))), hpv_sites$SV.id))
  
  # Make a matrix 
  
  ## TIME-CONSUMING STEP ##
  for (i in rownames(mat)) {
    # Find the events that the read spans and mark them as TRUE in the matrix
    events <- c(names(hpv_reads_filt[grep(i, hpv_reads_filt)]))
    mat[rownames(mat) == i, colnames(mat) %in% events] <- TRUE
  }
  
  mat[is.na(mat)] <- FALSE # label NAs as FALSE
  mat <- apply(mat, 2, as.logical)
  rownames(mat) <- unique(as.character(unlist(hpv_reads_filt)))
  
  # indicate which SV(s) the reads contain
  groups <- apply(mat, 1, which)
  groups <- lapply(groups, names)
  groups <- lapply(groups, function(x){paste(x, collapse= ",")})
  ugroups <- unique(unlist(groups)) # get the unique groups
  
  # make a conditional in the case that there are no overlapping sites
  if (all(ugroups == "")) {
    event_list <- as.list(hpv_sites$SV.id)
    names(event_list) <- paste0("event", 1:length(event_list)) 
    
  } else {
    ####
    #### COLLAPSE THE SITES INTO INTEGRATION EVENTS
    ####
    
    ugroups <- strsplit(x = ugroups, split = ",")
    names(ugroups) <- paste0("pattern", 1:length(ugroups))
    
    # Three rounds of collapsing the unique groups 
    mat2 <- matrix(nrow = length(ugroups), ncol = length(ugroups), dimnames = list(names(ugroups), names(ugroups)))
    mat3 <- matrix(nrow = length(ugroups), ncol = length(ugroups), dimnames = list(names(ugroups), names(ugroups)))
    mat4 <- matrix(nrow = length(ugroups), ncol = length(ugroups), dimnames = list(names(ugroups), names(ugroups)))
    
    # Round 1
    ugroups2 <- list()
    for (p in names(ugroups)) {
      pat <- ugroups[[p]]
      
      compare <- unlist(lapply(ugroups, function(x){test <- any(pat %in% x)}))
      
      mat2[p,] <- compare
      
      condensed <- unique(c(unlist(ugroups[which(mat2[p,])])))
      
      ugroups2[[p]] <- condensed
    }
    
    # Round 2
    ugroups3 <- list()
    for (p in names(ugroups2)) {
      pat <- ugroups2[[p]]
      
      compare <- unlist(lapply(ugroups2, function(x){test <- any(pat %in% x)}))
      
      mat3[p,] <- compare
      
      condensed <- unique(c(unlist(ugroups2[which(mat3[p,])])))
      
      ugroups3[[p]] <- condensed
    }
    
    # Round 3
    ugroups4 <- list()
    for (p in names(ugroups3)) {
      pat <- ugroups3[[p]]
      
      compare <- unlist(lapply(ugroups3, function(x){test <- any(pat %in% x)}))
      
      mat4[p,] <- compare
      
      condensed <- unique(c(unlist(ugroups3[which(mat4[p,])])))
      
      ugroups4[[p]] <- condensed
    }
    
    # Shorten the list to only the unique groups
    events <- lapply(ugroups4, function(x){paste(x, collapse= ",")})
    events <- unique(unlist(events))
    event_list <- strsplit(x = events, split = ",")
    names(event_list) <- paste0("event", 1:length(event_list)) 
  }
  
  # Reads in the whole event
  event_reads <- NULL
  for (event in names(event_list)) {
    sub <- hpv_reads_filt[event_list[[event]]]
    event_reads[[event]] <- as.character(unlist(sub))
    #write(event_reads[[event]], file = paste0(outdir,"/", event, ".txt"), sep = "\t")
  }
  
  # put the event reads into a df
  event_reads_df <- lapply(event_reads, function(x){
    df <- data.frame(reads = x)
  })
  event_reads_df <- bind_rows(event_reads_df,.id = "event")
  
  ####
  #### GROUP DIFFERENT SITES TOGETHER BY DISTANCE (500kb)
  ####
  bed <- data.frame(chr=hpv_sites$V1, start=hpv_sites$V2, end=hpv_sites$V2+1, site=hpv_sites$SV.id)
  bed_sort <- bedtoolsr::bt.sort(bed)
  b.merge <- bedtoolsr::bt.merge(i = bed_sort, d = 500000, o = "collapse", c = 4)
  
  dist_events <- strsplit(b.merge$V4, ",")
  
  ####
  #### COMBINE EVENTS THAT ARE PAIRED BY READ OR BY DISTANCE
  ####
  
  combo_mat <- matrix(nrow = length(hpv_sites$SV.id), ncol = length(hpv_sites$SV.id), dimnames = list(hpv_sites$SV.id, hpv_sites$SV.id))
  
  for (site in hpv_sites$SV.id){
    
    which.r <- lapply(event_list, function(x){is.na(which(x == site))})
    event.r <- grep(FALSE, which.r, fixed = T)
    
    which.d <- lapply(dist_events, function(x){is.na(which(x == site))})
    event.d <- grep(FALSE, which.d, fixed = T)
    
    e1 <- event_list[[event.r]]
    e2 <- dist_events[[event.d]]
    
    # combine together
    e <- unique(c(e1, e2))
    
    # fill the matrix
    combo_mat[site, e] <- TRUE
    combo_mat[site, hpv_sites$SV.id[!hpv_sites$SV.id %in% e]] <- FALSE
    
  }
  
  # indicate which SV(s) the reads contain
  final_groups <- apply(combo_mat, 1, function(x){list(names(which(x)))})
  final_groups <- lapply(final_groups, unlist)
  final_groups <- lapply(final_groups, function(x){paste(x, collapse= ",")})
  final_ugroups <- unique(unlist(final_groups)) # get the unique groups
  
  final_events1 <- strsplit(final_ugroups, ",")
  
  names(final_events1) <- paste0("pattern", 1:length(final_events1))
  
  # collapsing the unique groups 
  m <- matrix(nrow = length(final_events1), ncol = length(final_events1), dimnames = list(names(final_events1), names(final_events1)))
  
  final_events2 <- list()
  for (p in names(final_events1)) {
    pat <- final_events1[[p]]
    
    compare <- unlist(lapply(final_events1, function(x){test <- any(pat %in% x)}))
    
    m[p,] <- compare
    
    condensed <- unique(c(unlist(final_events1[which(m[p,])])))
    
    final_events2[[p]] <- condensed
  }
  
  # Shorten the list to only the unique groups
  final_events <- lapply(final_events2, function(x){paste(x, collapse= ",")})
  final_events <- unique(unlist(final_events))
  final_events <- strsplit(x = final_events, split = ",")
  names(final_events) <- paste0("event", 1:length(final_events)) 
  
  print(paste0("This sample has ", length(final_events), " HPV integration event(s)"))
  
  ####
  #### MAKE A DATAFRAME TO SAVE WITH SITE/EVENT INFO
  ####
  
  events <- data.frame(chr=NA, pos=NA, HPVchr=NA, HPVpos=NA, HPV.site=hpv_sites$SV.id, event=NA, VAF=NA)
  
  suppressMessages(for (i in 1:nrow(events)) {
    site <- events$HPV.site[i]
    events$event[i] <- names(final_events)[grep(paste0("\\b", site, "\\b"), final_events)]
  })
  
  events$chr <- hpv_sites$V1[match(events$HPV.site, hpv_sites$SV.id)]
  events$pos <- hpv_sites$V2[match(events$HPV.site, hpv_sites$SV.id)]
  events$HPVchr <- hpv_sites$chr2[match(events$HPV.site, hpv_sites$SV.id)]
  events$HPVpos <- hpv_sites$end[match(events$HPV.site, hpv_sites$SV.id)]
  events$VAF <- hpv_sites$VAF[match(events$HPV.site, hpv_sites$SV.id)]
  rownames(events) <- 1:nrow(events)
  
  events <- arrange(events, event, HPV.site)
  
  ####
  #### WRITE AND SAVE READ NAMES AND DATAFRAME
  ####
  
  # new event reads
  event_reads2 <- NULL
  for (event in names(final_events)) {
    sub <- hpv_reads_filt[final_events[[event]]]
    event_reads2[[event]] <- as.character(unlist(sub))
    write(event_reads2[[event]], file = paste0(outdir,"/", event, ".txt"), sep = "\t")
  }
  
  # Reads in individual sites
  for (site in names(hpv_reads_filt)) {
    write(hpv_reads_filt[[site]], file = paste0(outdir,"/", site, ".txt"), sep = "\t")
  }
  
  # add the number of reads to the events data frame
  num_reads <- stack(hpv_reads_filt)
  num_reads <- num_reads %>%
    group_by(ind) %>%
    summarise(n = n())
  
  events$read.depth <- num_reads$n[match(events$HPV.site, num_reads$ind)]
  
  # save the events table
  write.table(events, file = paste0(outdir,"/summary.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
  
  ####
  #### Write a bed file for the integration sites
  ####
  events_bed <- events[,c(1,2,2,5)]
  events_bed$pos.1 <- events_bed$pos.1 + 1
  write.table(events_bed, file = paste0(outdir,"/hpv_integration_sites.bed"), sep = "\t", quote = F, col.names = F, row.names = F)
  
}
