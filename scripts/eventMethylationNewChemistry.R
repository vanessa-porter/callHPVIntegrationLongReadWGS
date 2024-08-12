#!/gsc/software/linux-x86_64-centos7/R-4.0.2/bin/Rscript --vanilla
.libPaths("/projects/vporter_prj/R/x86_64-centos7-linux-gnu-library/4.0")

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
#me <- read.delim("/projects/hpv_nanopore_prj/htmcp/call_integration/output/HTMCP-03-06-02219/methylation/hpv_reads_methylation.tsv", comment.char = '#', header = T, sep = "\t", stringsAsFactors = F)
me <- read.delim(opt$methyl, comment.char = '#', header = T, sep = "\t", stringsAsFactors = F)

#dir <- "/projects/hpv_nanopore_prj/htmcp/call_integration/output/HTMCP-03-06-02219/events"
dir <- opt$dir

# event reads
files <- grep("event", list.files(dir, recursive = T, full.names = F),value = T)
files <- files[grep(".txt", files)]
files <- files[!grepl("dist", files)]
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

me$event <- reads$event[match(me$read_id, reads$read.name)]
me <- me[complete.cases(me),]
me <- me[me$ref_position > 0,]
me$methylation.call <- ifelse(me$call_code == "m", "methylated", "unmethylated")

freq <- me %>%
  group_by(event, chrom, ref_position, methylation.call) %>%
  summarise(n.calls = n()) %>%
  ungroup() %>%
  spread(methylation.call,  n.calls) %>%
  mutate(methylated = replace_na(methylated, 0),
         unmethylated = replace_na(unmethylated, 0))

freq$total.calls <- freq$methylated + freq$unmethylated
freq$perc.methylated <- freq$methylated / freq$total.calls
freq <- freq[freq$total.calls > 2,]

### ----------------------------------------------------------------
### Save
### ----------------------------------------------------------------

freq <- freq[,c(2,3,1,4:7)] %>% arrange(event, chrom, ref_position)
colnames(freq)[c(1,2)] <- c("chromosome", "start")
write.table(freq, file = opt$out, quote = F, col.names = T, row.names = F, sep = "\t")
#write.table(freq, file = "/projects/hpv_nanopore_prj/htmcp/call_integration/F46073/methylation/event_methyl_freq.tsv", quote = F, col.names = T, row.names = F, sep = "\t")
