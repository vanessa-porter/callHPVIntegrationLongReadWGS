#!/gsc/software/linux-x86_64-centos7/R-4.0.2/bin/Rscript --vanilla

#Note: these packages need to be installed.
suppressMessages(library(optparse))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))

# Make help options
option_list = list(
  make_option(c("-m", "--methyl"), type="character", default=NULL,
              help="HPV-selected methylation", metavar="character"),
  make_option(c("-d", "--dir"), type="character", default=NULL,
              help="HPV event directory", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="event_methylation.tsv",
              help="Output file name", metavar="character")
)

# load in required 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

### ----------------------------------------------------------------
### Import the files
### ----------------------------------------------------------------

# HPV reads methylation file
me <- read.delim(opt$methyl, comment.char = '#', header = F, sep = "\t", stringsAsFactors = F)
dir <- opt$dir

# event reads
files <- grep("event", list.files(dir, recursive = T, full.names = F),value = T)
files <- files[grep(".txt", files)]
files <- files[!grepl("depth", files)]
files <- grep(paste(files,collapse="|"), list.files(dir, recursive = T, full.names = T),value = T)
names <- gsub(paste0(dir, "/"), "", files)
names <- gsub(".txt", "", names)

event_reads <- NULL
for (i in 1:length(files)) {
  event_reads[[i]] <- read.delim(files[i], header = F, sep = "\t", stringsAsFactors = F)
}
names(event_reads) <- names
reads <- bind_rows(event_reads, .id = "event")
colnames(reads)[2] <- "read.name"

### ----------------------------------------------------------------
### Separate methylation frequencies by event
### ----------------------------------------------------------------

colnames(me) <- c("chromosome","strand","start","end","read_name","log_lik_ratio",
                  "log_lik_methylated","log_lik_unmethylated","num_calling_strands","num_motifs","sequence")

me$event <- reads$event[match(me$read_name, reads$read.name)]
me <- me[complete.cases(me),]
me$methylation.call <- ifelse(me$log_lik_ratio > 0, "methylated", "unmethylated")

freq <- me %>%
  group_by(event, chromosome, start, methylation.call) %>%
  summarise(n.calls = n()) %>%
  ungroup() %>%
  spread(methylation.call,  n.calls) %>%
  mutate(methylated = replace_na(methylated, 0),
         unmethylated = replace_na(unmethylated, 0))

freq$total.calls <- freq$methylated + freq$unmethylated
freq$perc.methylated <- freq$methylated / freq$total.calls
freq <- freq[freq$total.calls > 4,]

### ----------------------------------------------------------------
### Save
### ----------------------------------------------------------------

freq <- freq[,c(2,3,1,4:7)] %>% arrange(event, chromosome, start)
write.table(freq, file = opt$out, quote = F, col.names = T, row.names = F, sep = "\t")
