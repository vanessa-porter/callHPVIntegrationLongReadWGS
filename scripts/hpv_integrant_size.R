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


### -------------------------------------------------------------------------------
### READ IN FILES - best example = HTMCP-03-06-02238
### -------------------------------------------------------------------------------

#pafr <- read_paf("/projects/hpv_nanopore_prj/htmcp/call_integration/output/HTMCP-03-06-02097/bam/hpv_reads.paf")
#sum <- read.delim("/projects/hpv_nanopore_prj/htmcp/call_integration/output/HTMCP-03-06-02097/events/summary.txt", header = T)
#out <- "/projects/hpv_nanopore_prj/htmcp/call_integration/output/HTMCP-03-06-02097/hpv_size/"

hpv <- read.delim("/projects/hpv_nanopore_prj/refs/HPV_all_chromsizes.txt", header = F)

# Make help options
option_list = list(
  make_option(c("-p", "--paf"), type="character", default=NULL,
              help="paf file with HPV reads", metavar="character"),
  make_option(c("-s", "--sum"), type="character", default=NULL,
              help="summary file for the integration events", metavar="character"),
  make_option(c("-o", "--out"), type="character", default = "./",
              help="Output folder name", metavar="character")
)

# load in options 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
out <- opt$out
pafr <- read_paf(opt$paf)
sum <- read.delim(opt$sum, header = T)

### -------------------------------------------------------------------------------
### READ IN FILES - best example = HTMCP-03-06-02238
### -------------------------------------------------------------------------------

# ggplot directly
# Flip coordinates for (-) strand:
paf.tab <- pafr %>%
  mutate(temp = tstart) %>%
  mutate(tstart=ifelse(strand=='-', tend, tstart)) %>%
  mutate(tend=ifelse(strand=='-', temp, tend)) %>%
  mutate(temp = tstart) %>%
  filter(mapq == 60) %>%
  #mutate(tstart=ifelse(strand=='-', tend, tstart)) %>%
  #mutate(tend=ifelse(strand=='-', temp, tend)) %>%
  #mutate(qname = fct_reorder(qname, desc(qlen))) %>%
  dplyr::select(1:11) %>%
  arrange(qname, qstart)

paf.tab <- paf.tab[!duplicated(paf.tab),]
rownames(paf.tab) <- 1:nrow(paf.tab)

### -------------------------------------------------------------------------------
### REOCCURING BREAKS
### -------------------------------------------------------------------------------

# make a range of values that each break could be (+/- 5bp) 
# the calls are imperfect so some leniency is added so we can get more data out of the reads
sumBrks <- sum$pos

brk_range <- list()
for (i in 1:length(sumBrks)) {
  n <- sumBrks[i]
  range <- c((n-30):(n+30))
  brk_range[[i]] <- data.frame(possible_values = as.numeric(range))
}
names(brk_range) <- sumBrks
brk_range <- bind_rows(brk_range, .id = "breakpoint")

# make an alternate tstart/tend with the matching breakpoints (human and HPV)
paf.tab$tstart2 <- NA
paf.tab$tend2 <- NA

for (i in 1:nrow(paf.tab)) {
  sval <- as.numeric(paf.tab[i,"tstart"])
  eval <- as.numeric( paf.tab[i,"tend"])
  sbr <- brk_range$breakpoint[brk_range$possible_values == sval]
  ebr <- brk_range$breakpoint[brk_range$possible_values == eval]
  if(length(sbr) > 0){
    paf.tab$tstart2[i] <- brk_range$breakpoint[brk_range$possible_values == sval]
  }
  if(length(ebr) > 0){
    paf.tab$tend2[i] <- brk_range$breakpoint[brk_range$possible_values == eval]
  }
}

# make a range of values that each break could be (+/- 5bp) -- HPV BREAKS
sumHPVBrks <- sum$HPVpos

brk_hpv_range <- list()
for (i in 1:length(sumHPVBrks)) {
  n <- sumHPVBrks[i]
  range <- c((n-30):(n+30))
  brk_hpv_range[[i]] <- data.frame(possible_values = as.numeric(range))
}
names(brk_hpv_range) <- sumHPVBrks
brk_hpv_range <- bind_rows(brk_hpv_range, .id = "breakpoint")

# make an alternate tstart/tend with the matching breakpoints (human and HPV)
brk <- rbind.data.frame(brk_hpv_range, brk_range)
paf.tab$tstart2 <- NA
paf.tab$tend2 <- NA

for (i in 1:nrow(paf.tab)) {
  sval <- as.numeric(paf.tab[i,"tstart"])
  eval <- as.numeric( paf.tab[i,"tend"])
  sbr <- brk_range$breakpoint[brk$possible_values == sval]
  ebr <- brk_range$breakpoint[brk$possible_values == eval]
  if(length(sbr) > 0){
    paf.tab$tstart2[i] <- brk$breakpoint[brk$possible_values == sval]
  }
  if(length(ebr) > 0){
    paf.tab$tend2[i] <- brk$breakpoint[brk$possible_values == eval]
  }
}

print("check1")

### -------------------------------------------------------------------------------
### GET THE LENGTH OF COMPLETE HPV INTEGRATIONS
### GET THE MAXIMUM LENGTH OF INCOMPLETE HPV INTEGRATIONS
### -------------------------------------------------------------------------------

reads <- unique(paf.tab$qname)
rList <- list()
for (r in 1:length(reads)) {
  p <- as.data.frame(paf.tab[paf.tab$qname == reads[r],])
  p$breakpoint <- ifelse(grepl("chr",p$tname), as.numeric(rownames(p)), "HPV") 
  p$breakpoint2 <- NA
  chrFactors <- as.factor(p$breakpoint[p$breakpoint != "HPV"])
  chrFactors <- factor(chrFactors, levels = chrFactors)
  p$breakpoint2[p$breakpoint != "HPV"] <- paste0("break",as.numeric(chrFactors))
  breaks <- unique(p$breakpoint2[!is.na(p$breakpoint2)])
  order <- as.numeric(rownames(p[p$breakpoint2 %in% breaks,]))
  
  if(!all(is.na(p$tstart2)) & !all(is.na(p$tend2))){
      if(length(breaks) > 1 & sum(!is.na(p$tend2)) > 1 & sum(!is.na(p$tstart2))){
          df <- data.frame(read = rep(reads[r],length(2:length(breaks))), 
                           break1 = NA, break2 = NA, 
                           HPVbreak1 = NA, HPVbreak2 = NA,  
                           length = NA,
                           status=rep("complete",length(2:length(breaks))))
          for (i in 2:length(breaks)) {
              break1 <- breaks[i-1]
              break2 <- breaks[i]
              
              # get the rows in between break1 and break2
              psub1 <- p[p$breakpoint2 == break1,]
              psub1 <- psub1[complete.cases(psub1[,-c(12,13)]),]
              rows1 <- as.numeric(rownames(psub1))
              psub2 <- p[p$breakpoint2 == break2,]
              psub2 <- psub2[complete.cases(psub2[,-c(12,13)]),]
              rows2 <- as.numeric(rownames(psub2))
              btw <- c((rows1+1):(rows2-1))
              
              # find which chromosome is between the breaks
              btw_chr <- unique(p$tname[rownames(p) %in% btw])
              
              if(grepl("HPV", btw_chr)){
                  pbtw <- p[rownames(p) %in% btw,]
                  btw_start <- p$tstart2[rownames(p) == min(btw)]
                  btw_end <- p$tend2[rownames(p) == max(btw)]
                  length <-  max(pbtw$qend) - min(pbtw$qstart)
                  
                  # make a dataframe with the length of HPV in between the breaks
                  if(!is.na(psub1$tend2) & !is.na(psub2$tstart2) & !is.na(btw_start) & !is.na(btw_end) & psub2$tstart2 != psub1$tend2){
                      bps <- c(as.numeric(psub1$tend2), as.numeric(psub2$tstart2))
                      chrs <- c(psub1$tname, psub2$tname)
                      hpv_bps <- c(btw_start,btw_end)
                      
                      b1 <- which.min(bps)
                      b2 <- which.max(bps)
                      
                      df[i-1,] <- data.frame(read = reads[r], 
                                             break1 = paste0(chrs[b1], ":", bps[b1]), 
                                             break2 = paste0(chrs[b2], ":", bps[b2]), 
                                             HPVbreak1 = paste0(btw_chr, ":", hpv_bps[b1]),  
                                             HPVbreak2 = paste0(btw_chr, ":", hpv_bps[b2]),  
                                             length = length, status = "complete")
                      df <- df[complete.cases(df),]
                      rList[[r]] <- df
                  }
              }
          }
      } else if (length(breaks) == 1 & nrow(p) > 1){
          # find the length of the HPV chromosomes following the break
          bp <- rownames(p)[p$breakpoint2 == breaks & complete.cases(p$breakpoint2)]
          other_rows <- rownames(p)[is.na(p$breakpoint2)]
          
          # test that the hpv alignments before/after the break are sequential 
          diff_test <- rle(diff(as.numeric(other_rows)))
          
          # get the positions of the HPV alignments with the max * of sequential alignments
          if(length(diff_test$lengths) == 3){
              if(diff_test$lengths[1] > diff_test$lengths[3]){
                  other_rows <- other_rows[1:(diff_test$lengths[1]+1)]
              } else {
                  other_rows <- other_rows[!1:length(other_rows) %in% 1:(diff_test$lengths[1]+1)]
              }
          } else if(length(diff_test$lengths) == 2){
              run <- which(diff_test$values == 1)
              
              if(run ==1){
                  other_rows <- other_rows[1:(diff_test$lengths[1]+1)]
              } else if(run ==2){
                  other_rows <- other_rows[!1:length(other_rows) %in% 1:(diff_test$lengths[1])]
              }
          }
          
          # test if the position of the HPV break has a reoccurring breakpoint
          hpv_break <- other_rows[which((abs(as.numeric(bp) - as.numeric(other_rows))) == 1)]
          hpv_start <- p[hpv_break,"tstart2"]
          hpv_end <- p[hpv_break,"tend2"]
          
          # get the hpv breakpoint depending on if the alignments are before or after the human alignment
          if(as.numeric(bp) > as.numeric(hpv_break)){
              hpvbreak <- hpv_end
          } else if(as.numeric(bp) < as.numeric(hpv_break)){
              hpvbreak <- hpv_start
          }
          
          # get the breakpoint of the human chromosome
          bp_start <- p[bp,"tstart2"]
          bp_end <- p[bp,"tend2"]
          
          if(as.numeric(bp) > as.numeric(hpv_break)){
              bp_pos <- bp_start
          } else if(as.numeric(bp) < as.numeric(hpv_break)){
              bp_pos <- bp_end
          }
          
          # get the length of the incomplete HPV segment
          phpv <- p[rownames(p) %in% other_rows,]
          length <-  max(phpv$qend) - min(phpv$qstart)
          
          # get the breakpoint position
          break1 = paste0(p[bp,"tname"], ":", bp_pos)
          
          # put the values in the dataframe
          if(!is.na(hpvbreak)){
              hpv_pos <- paste0(p[hpv_break, "tname"],":",hpvbreak)
              df <- data.frame(read = reads[r], 
                               break1 = break1, break2 = NA, HPVbreak1 = hpv_pos, HPVbreak2 = NA,  length = length, status="incomplete")
          }
          
          rList[[r]] <- df
          
      } else if(length(breaks) > 1 & max(diff(order)) == 1){
          # find the length of the HPV chromosomes following the breaks
          bp <- rownames(p)[p$breakpoint2 == breaks & complete.cases(p$breakpoint2)]
          other_rows <- rownames(p)[is.na(p$breakpoint2)]
          
          # test that the hpv alignments before/after the break are sequential 
          diff_test <- rle(diff(as.numeric(other_rows)))
          
          # get the positions of the HPV alignments with the max * of sequential alignments
          if(length(diff_test$lengths) == 3){
              if(diff_test$lengths[1] > diff_test$lengths[3]){
                  other_rows <- other_rows[1:(diff_test$lengths[1]+1)]
              } else {
                  other_rows <- other_rows[!1:length(other_rows) %in% 1:(diff_test$lengths[1]+1)]
              }
          } else if(length(diff_test$lengths) == 2){
              run <- which(diff_test$values == 1)
              
              if(run ==1){
                  other_rows <- other_rows[1:(diff_test$lengths[1]+1)]
              } else if(run ==2){
                  other_rows <- other_rows[!1:length(other_rows) %in% 1:(diff_test$lengths[1])]
              }
          }
          
          # test if the position of the HPV break has a reoccurring breakpoint
          hpv_break <- other_rows
          hpv_start <- p[other_rows,"tstart2"]
          hpv_end <- p[other_rows,"tend2"]
          
          # get the hpv breakpoint depending on if the alignments are before or after the human alignment
          if(min(as.numeric(bp)) > max(as.numeric(other_rows))){
              hpvbreak <- hpv_end
          } else if(max(as.numeric(bp)) < min(as.numeric(other_rows))){
              hpvbreak <- hpv_start
          }
          
          # get the breakpoint of the human chromosome
          bp_start <- p[bp,"tstart2"]
          bp_end <- p[bp,"tend2"]
          
          if(min(as.numeric(bp)) > max(as.numeric(other_rows))){
              bp_pos <- bp_start
          } else if(max(as.numeric(bp)) < min(as.numeric(other_rows))){
              bp_pos <- bp_end
          }
          
          # get the length of the incomplete HPV segment
          phpv <- p[rownames(p) %in% other_rows,]
          length <-  max(phpv$qend) - min(phpv$qstart)
          
          # get the breakpoint position
          break1 = paste0(p[bp,"tname"], ":", bp_pos)
          
          # put the values in the dataframe
          if(!all(is.na(hpvbreak))){
              hpv_pos <- paste0(p[hpv_break, "tname"],":",hpvbreak)
              df <- data.frame(read = reads[r], 
                               break1 = break1, break2 = NA, HPVbreak1 = hpv_pos, HPVbreak2 = NA,  length = length, status="incomplete")
          }
          
          rList[[r]] <- df
    }
  }
}

mylist <- rList[lengths(rList) != 0]
mylist <- lapply(mylist, as.data.frame)
hpvSize <- bind_rows(mylist)
hpvSize <- hpvSize[complete.cases(hpvSize$break1),]

hpvSize <- hpvSize[!grepl("NA", hpvSize$break1),]
hpvSize <- hpvSize[!grepl("NA", hpvSize$HPVbreak1),]
hpvSize$bp_pair <- paste0(hpvSize$break1, "_", hpvSize$break2)

print("check2")

### -------------------------------------------------------------------------------
### GET THE SIZES AND SUMMARY OF EACH BREAKPOINT PAIR
### -------------------------------------------------------------------------------

hpvlength <- hpv$V2[hpv$V1 == unique(sum$HPVchr)]
hpvSize$n_HPV <- hpvSize$length / hpvlength
hpvSize$bp_pair_name <- paste0("break-pair",as.numeric(as.factor(hpvSize$bp_pair)))
hpvSize$event <- sum$event[match(hpvSize$break1, paste0(sum$chr, ":", sum$pos))]
hpvSizeComp <- hpvSize[hpvSize$status == "complete",]
hpvSizeIncomp <- hpvSize[hpvSize$status == "incomplete",]

### -------------------------------------------------------------------------------
### GROUP THE LENGTHS PER BP PAIR
### -------------------------------------------------------------------------------

bp_pair_names <- unique(hpvSizeComp$bp_pair_name)

if(length(bp_pair_names) > 0){
    lgroups <- list()
    for(j in 1:length(bp_pair_names)){
        # subset the dataframe
        bp_pair_name <- unique(hpvSizeComp$bp_pair_name)[j]
        sub <- hpvSizeComp[hpvSizeComp$bp_pair_name == bp_pair_name,]
        sub <- arrange(sub, length)
        rownames(sub) <- 1:nrow(sub)
        
        # how many groups of lengths
        d <- diff(sub$length) > 300
        n_lgroups <- sum(d) + 1
        
        if (n_lgroups == 1){
            gl <- sub$length
            g <- "group1"
            lg <- data.frame(group = rep(g, length(gl)), length = gl)
        } else if (n_lgroups == 2){
            pos <- which(d == T)
            l <- sub$length
            
            ## first group
            start1 <- 1 
            end1 <- pos[1]
            gl1 <- l[start1:end1]
            g1 <- c("group1")
            lg1 <- data.frame(group = rep(g1, length(gl1)), length = gl1)
            
            ## second group
            start2 <- pos + 1
            end2 <- length(l)
            gl2 <- l[start2:end2]
            g2 <- c("group2")
            lg2 <- data.frame(group = rep(g2, length(gl2)), length = gl2)
            
            lg <- rbind.data.frame(lg1, lg2)
        }else {
            pos <- which(d == T)
            lg <- list()
            for (i in 1:n_lgroups) {
                l <- sub$length
                if(i == 1){ # if the first group
                    if(d[1] == T){
                        gl <- l[1] 
                    } else {
                        start <- 1 
                        end <- pos[1]
                        gl <- l[start:end]
                    }
                } else if (i == n_lgroups){ # if the last group
                    if(d[length(d)] == T){
                        gl <- l[length(l)]
                    } else{
                        start <- pos[length(pos)]
                        end <- length(l)
                        gl <- l[start:end]
                    }
                } else { # all other groups
                    start <- pos[i-1] + 1
                    end <- pos[i]
                    gl <- l[start:end]
                }
                
                # put dataframe together
                g <- paste0("group", i)
                df <- data.frame(group = rep(g, length(gl)), length = gl)
                lg[[i]] <- df 
            }
            lg <- bind_rows(lg)
        }
        lg$bp_pair_name <- bp_pair_name
        lg$bp_pair <- unique(sub$bp_pair)
        lg <- lg[,c(4,3,1,2)]
        lgroups[[j]] <- lg
    }
    
    lgroups <- bind_rows(lgroups)
    lgroups$n_HPV <- lgroups$length / hpvlength
    lgroups$bp_pair_group <- paste0(lgroups$bp_pair_name, "_", lgroups$group)
    lgroups$status <- "complete"
    lgroups$event <- hpvSize$event[match(lgroups$bp_pair_name, hpvSize$bp_pair_name)]
    
    ### -------------------------------------------------------------------------------
    ### SUMMARIZE THE LENGTH GROUPS
    ### -------------------------------------------------------------------------------
    
    sizeSum <- lgroups %>%
        group_by(bp_pair_group) %>%
        summarise(length = median(length),
                  nHPV = median(n_HPV),
                  nreads = n(),
                  bp_pair = unique(bp_pair),
                  bp_pair_name = unique(bp_pair_name)) %>%
        mutate(status = "complete") %>%
        select(bp_pair_name, bp_pair, bp_pair_group, length, nHPV, nreads, status)
    
    sizeSum2 <- sizeSum %>%
        group_by(bp_pair_name) %>%
        summarise(diff = paste(round(diff(nHPV),digits = 1),sep = ","))
    
    sizeSum3 <- sizeSum %>%
        group_by(bp_pair_name) %>%
        summarise(ngroups = n())
    
    sizeSum$ngroups <- sizeSum3$ngroups[match(sizeSum$bp_pair_name, sizeSum3$bp_pair_name)]
    
    ### -------------------------------------------------------------------------------
    ### FIGURES
    ### -------------------------------------------------------------------------------
    
    plot <- ggplot(lgroups, aes(x = bp_pair, y = n_HPV, colour = group)) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 3) +
        scale_y_continuous(
            # Features of the first axis
            name = "# of HPV Genomes",
            # Add a second axis and specify its features
            sec.axis = sec_axis( trans=~.*hpvlength, name="Length of HPV Integrant"),
            expand = c(0, 0), limits = c(0, NA),
        ) +
        theme_bw() +
        scale_colour_lancet() +
        facet_grid(~ event, scales = "free_x", space = "free")+
        labs(x = "Breakpoint Pair") +
        theme(axis.title = element_text(size = 14, colour = "black", face = "bold"),
              axis.text.y = element_text(size = 12, colour = "black"),
              axis.text.x = element_text(size = 10, colour = "black", angle = 60, hjust = 1, vjust = 1),
              panel.grid.minor = element_blank(),
              strip.text = element_text(size = 12, colour = "black"),
              legend.text = element_text(size = 12, colour = "black"),
              legend.title = element_text(size = 14, colour = "black", face = "bold"))
    
    ### -------------------------------------------------------------------------------
    ### SAVE TABLES AND FIGURES TO THE SPECIFIED DIRECTORY
    ### -------------------------------------------------------------------------------
    
    write.table(lgroups, file = paste0(out, "/hpvLengthGroups.txt"), quote = F, col.names = T, row.names = F, sep = "\t")
    write.table(sizeSum, file = paste0(out, "/hpvSizeSummary.txt"), quote = F, col.names = T, row.names = F, sep = "\t")
    ggsave(plot, filename = paste0(out, "/hpvIntegrantSize.pdf"), width = 10, height = 6, units = "in")
}

### -------------------------------------------------------------------------------
### INCOMPLETE HPV INTEGRANTS
### -------------------------------------------------------------------------------

sizeIncomSum <- hpvSizeIncomp %>%
  filter(!break1 %in% c(unique(hpvSizeComp$break1),unique(hpvSizeComp$break2))) %>%
  group_by(break1, bp_pair_name) %>%
  summarise(length = max(length),
            nHPV = max(n_HPV),
            nreads = n()) %>%
  mutate(group = NA, 
         bp_pair_group = break1,
         status = "incomplete",
         ngroups = 0) %>%
  select(bp_pair_name, break1, bp_pair_group, length, nHPV, nreads, status, ngroups)

colnames(sizeIncomSum)[2] <- "bp_pair"

### -------------------------------------------------------------------------------
### CATEGORIZE THE BP PAIRS 
### CATEGORIES = PARTIAL, FULL, HETEROLOGOUS, INCOMPLETE
### HPV SIZE = PARTIAL, 1+, 2+, 3+ (BASED ON THE HIGHEST NUMBER IN HETEROLOGOUS BP PAIRS)
### -------------------------------------------------------------------------------

# combine dfs
if (length(bp_pair_names) > 0){
    size <- rbind.data.frame(sizeSum, sizeIncomSum)
} else{
    size <- rbind.data.frame(sizeIncomSum)
}

categories <- size %>%
  group_by(bp_pair_name,bp_pair) %>%
  summarise(ngroups = unique(ngroups),
            max_length = max(length),
            max_nHPV = max(nHPV),
            status = unique(status)) %>%
  mutate(
    category = case_when(
      status == "incomplete" ~ "incomplete",
      ngroups > 1 ~ "heterologous",
      max_nHPV > 1 ~ "full",
      TRUE ~ "partial"),
    size_category = case_when(
      max_nHPV > 3 ~ "over3",
      max_nHPV > 2 ~ "over2",
      max_nHPV > 1 ~ "over1",
      TRUE ~ "less1"))


### -------------------------------------------------------------------------------
### SAVE TABLES AND FIGURES TO THE SPECIFIED DIRECTORY
### -------------------------------------------------------------------------------

# save objects: hpvSize, lgroups, categories, sizeSum, plot

write.table(hpvSize, file = paste0(out, "/hpvSizeReads.txt"), quote = F, col.names = T, row.names = F, sep = "\t")
write.table(categories, file = paste0(out, "/hpvSizeCategories.txt"), quote = F, col.names = T, row.names = F, sep = "\t")

### -------------------------------------------------------------------------------
### DOT PLOT OF A SELECTED READ
### -------------------------------------------------------------------------------
#bp <- paf.tab[paf.tab$qname %in% hpvSize$read[hpvSize$bp_pair_name == "break-pair18"],]
  
### plot
#ggplot(bp %>% filter(tname %in% c(unique(sum$HPVchr), unique(sum$chr)))) +
#  aes(x=qstart, xend=qend, y=tstart, yend=tend, colour=tname) +
#  geom_segment(size = 1) +
#  xlab("Read Position") + ylab("Reference Position") +
#  theme_bw() +
#  scale_color_lancet() +
  #geom_hline(yintercept = 189748099, linetype = 2) +
#  facet_grid(tname ~ qname, shrink=TRUE, scales = "free", space = "free") + 
#  theme(panel.grid = element_blank())




