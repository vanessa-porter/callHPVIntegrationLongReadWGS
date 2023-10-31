#!/gsc/software/linux-x86_64-centos7/R-4.0.2/bin/Rscript --vanilla

#Note: these packages need to be installed.
library(optparse)

option_list = list(
  make_option(c("-1", "--methfile1"), type="character", default=NULL, 
              help="First methylation frequencies file", metavar="character"),
  make_option(c("-2", "--methfile2"), type="character", default=NULL, 
              help="Second methylation frequencies file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="Output file name 1", metavar="character"),
  make_option(c("-t", "--threads"), type="numeric", default=32, 
              help="Number of CPU threads to use [default= %default]", metavar="numeric"),
  make_option(c("-m", "--minreads"), type="numeric", default=3, 
              help="Minimum number of reads per CpG site to be included [default= %default]", metavar="numeric"),
  make_option(c("-f", "--fai"), type="character", 
              help="FASTA index for ref genome for splitting [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$methfile1)||is.null(opt$methfile2)||is.null(opt$out))
{
  stop("Both input files and the output file must be specified. Please see usage (--help).")
}

threads <- opt$threads
#threads <- 40

suppressMessages(library(DSS))
suppressMessages(library(data.table))

min_reads <- opt$minreads
#min_reads <- 3


#Load data and convert to DSS' format:

hap.1 <- fread(opt$methfile1)

hap_1 <- hap.1[,c(1,2,5,6)]
colnames(hap_1) <- c('chr', 'pos', 'N', 'X')
hap_1 <- hap_1[which(hap_1$N >= min_reads),]
hap_1 <- hap_1[hap_1$chr %in% c(paste0("chr", 1:22), "chrX"),]

hap.2 <- fread(opt$methfile2)

hap_2 <- hap.2[,c(1,2,5,6)]
colnames(hap_2) <- c('chr', 'pos', 'N', 'X')
hap_2 <- hap_2[which(hap_2$N >= min_reads),]
hap_2 <- hap_2[hap_2$chr %in% c(paste0("chr", 1:22), "chrX"),]

DSObject <- makeBSseqData(list(hap_1, hap_2), c('hap1','hap2'))

#Remove the raw data to save memory:
rm(hap.1, hap.2, hap_1, hap_2)
gc()

#Get chr lengths so we can divide the genome into tiles for parallelisation:

fai.info <- fread(opt$fai)
fai.info <- as.data.frame(fai.info)
rownames(fai.info) <- fai.info$V1
seqlengths(DSObject) <- fai.info[names(seqlengths(DSObject)),'V2']

#Next, use BiocParallel to scatter computation over the genomic chunks:
tiles_for_parallel <- tileGenome(seqlengths(DSObject), ntile=threads)
DSObjects_for_parallel <- lapply(tiles_for_parallel, function(x,y){subsetByOverlaps(y,x)}, DSObject)

calcOneTileDMRs <- function(DSObject)
{
  asm_dml <- DMLtest(DSObject, 'hap1', 'hap2', smoothing=TRUE, BPPARAM = MulticoreParam(workers=1) )
  asm_dmr <- callDMR(asm_dml, p.threshold = 0.01)
  return(asm_dmr)
}

mParam = MulticoreParam(workers=threads, progressbar=TRUE)

all <- bplapply(DSObjects_for_parallel, calcOneTileDMRs, BPPARAM = mParam)
all_df <- lapply(all, as.data.frame)
bind <- do.call("rbind", all_df)

fwrite(bind, opt$out, sep='\t')
