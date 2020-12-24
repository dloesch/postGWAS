#!/usr/bin/Rscript

library(ggplot2)
library(data.table)
library(ggrepel)
library(cowplot)
library(argparser)
library(tools)

#set argument parser
argp <- arg_parser("plot locuszoom-style plot with local LD data", hide.opts = TRUE)

#reguired arguments
argp <- add_argument(argp, "--sumstats", help="GWAS summary statistics", type = "character")
argp <- add_argument(argp, "--geno", help="genotype file", type = "character")

#optional arguments
argp <- add_argument(argp, "--window", help="kb window for LD", type = "numeric", default=250)
argp <- add_argument(argp, "--prefix", help="output prefix", type = "character", default="output")
argp <- add_argument(argp, "--snp", help="SNP; defaults to snp with lowest p-value", default = NA, type="character")
argp <- add_argument(argp, "--snpcol", help="column name for snp IDs", type = "character", default = NA)
argp <- add_argument(argp, "--pvalcol", help="column name for p-values",type = "character", default = NA)
argp <- add_argument(argp, "--chromcol", help="column name for chromosome", type = "character", default = NA)
argp <- add_argument(argp, "--poscol", help="column name for bp positions", type = "character", default = NA)
argp <- add_argument(argp, "--build", help="genome build of data", type = "character", default = "HG38")

#parse arguments
argsv <- parse_args(argp)


#creage gene plot
build <- argsv$build
if(build == "HG38" | build == "hg38"){
  genes <- read.csv("genes.pos.HG38.csv.gz", header=TRUE)
}

if(build == "HG19" | build == "hg19" | build == "hg37" | build == "HG37"){
  genes <- read.csv("genes.pos.HG19.csv.gz", header=TRUE)
}
#genes <- fread("genes_by_pos.csv.gz") #throws error

#summary stats and genotype data
geno.file <- argsv$geno
sumstats.file <- argsv$sumstats

#check input and exit
if(is.na(geno.file) | is.na(sumstats.file)){

  stop("Please provide GWAS summary stats and corresponding genotype file")
}

data <- fread(sumstats.file, header=TRUE, data.table = FALSE)

#column names
snpid <- argsv$snpcol
if(is.na(snpid)){
  snpid <- colnames(data)[colnames(data) %in% c("SNP", "snp", "Snp", "variant", "VARIANT", "ID")]
  if(length(snpid) != 1){
    stop("unable to infer SNP ID column name; please provide")
  }
}

pval <- argsv$pvalcol
if(is.na(pval)){
  pval <- colnames(data)[colnames(data) %in% c("P", "PVALUE", "pvalue", "p", "Score.pval", "PVAL", "pval")]
  if(length(pval) != 1){
    stop("unable to infer pval column name; please provide")
  }
}

chrom <- argsv$chromcol
if(is.na(chrom)){
  chrom <- colnames(data)[colnames(data) %in% c("C", "CHROM", "CHR", "c", "chrom", "chr", "chromosome", "CHROMOSOME")]
  if(length(chrom) != 1){
    stop("unable to infer chrom column name; please provide")
  }
}

pos <- argsv$poscol
if(is.na(pos)){
  pos <- colnames(data)[colnames(data) %in% c("pos", "POS", "BP", "bp", "position", "POSITION")]
  if(length(pos) != 1){
    stop("unable to infer BP column name; please provide")
  }
}

#range desired in kb
win <- argsv$window

snp <- argsv$snp

if( is.na(snp)){
  snp <- data[[snpid]][data[[pval]] == min(data[[pval]])][1]
}

chr <- data[[chrom]][data[[snpid]] == snp]
chr <- as.numeric(sub("chr", "", chr))
snp.pos <- as.numeric(data[[pos]][data[[snpid]] == snp])
start <- snp.pos - win*1000
end <- snp.pos + win*1000

genes <- genes[genes$chrom == chr & genes$end > start & genes$start < end,]
genes <- genes[!duplicated(genes),]
genes$start[genes$start < start] <- start
genes$end[genes$end > end] <- end

gene_panel <- ggplot(data = genes) + 
  geom_linerange(aes(x = gene_id, ymin = start, ymax = end)) +
  coord_flip() + ylab("") + 
  geom_text(aes(x = gene_id, y = start, label = gene_id), fontface = 2, alpha = I(0.7), hjust = "right", size= 2.5) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        strip.text.y = element_text(angle = 0),
        legend.position="bottom", 
        panel.grid.major.y = element_blank(),
        panel.background = element_blank()) +
  theme_void()


#ld data
#make sure data does not have an r2 column
data$R2 <- NULL

#get appropriate genotype file
isVCF <- file_ext(sub(".gz", "", geno.file)) == "vcf"
isBCF <- file_ext(geno.file) == "bcf"
isBED <- file_ext(geno.file) == "bed"


win <- win*2

tempfile <- paste0("temp", sample(1:1E06,1))
if(isBCF | isVCF){
  temp.vcf <- paste0(tempfile, ".vcf.gz")
  system2("bcftools", 
          paste("view -r", paste0("chr",chr, ":", start, "-", end), "-Oz -o",temp.vcf, geno.file))
  
  system2("plink",
          paste("--vcf", temp.vcf, "--r2 --ld-window 10000 --ld-window-kb", win, "--ld-snp", snp, "--ld-window-r2 0.0 --out", tempfile),
          stdout =FALSE)
}else if(isBED){
  plink.file <- sub(".bed", "", geno.file)
  system2("plink",
          paste("--bfile", plink.file, "--r2 --ld-window 10000 --ld-window-kb", win, "--ld-snp", snp, "--ld-window-r2 0.0 --out", tempfile),
          stdout =FALSE)
}else{
  stop("genotype file needs to be a vcf, bcf, or bed file")
}



ld <- fread(paste0(tempfile,".ld"), header=TRUE, stringsAsFactors = FALSE, data.table = FALSE)

system2("rm", paste0(tempfile,"*"))


#parse association statistics
data <- data[data[[chrom]] == chr,]
data <- data[data[[pos]] > start & data[[pos]] < end,]


#add ld value
r2 <- ld[c("SNP_B", "R2")]
data <- merge(data, r2, by.x=snpid, by.y="SNP_B", all.x=TRUE)

data$R2 <- as.numeric(as.character(data$R2))


#rename columns
dat2 <- data[c(pos, pval, "R2")]
colnames(dat2) <- c("POS", "PVALUE", "R2")
dat2$PVALUE <- -log10(dat2$PVALUE)
dat2 <- dat2[complete.cases(dat2$R2),]

#plot locuszoom
prefix <- argsv$prefix
outfile <- paste0(prefix, ".pdf")

pdf(outfile)
lz <- ggplot(dat2, aes(x=POS, y=PVALUE, color=R2)) + geom_point() +
  scale_color_gradient2(midpoint=0.5, low="darkblue", mid="green", high="red", na.value = "black") +
  ylim(0,max(dat2$PVALUE)+0.75)
lz <- lz  + theme_classic() + labs(title=paste0("Chromosome ", chr,":",snp.pos," locus"),x="Position", y="- log10(p-value)")
lz <- lz + geom_hline(yintercept=-log10(5E-08), linetype="dashed", color = "black")
lz <- lz + geom_text_repel(data = subset(dat2, PVALUE == max(dat2$PVALUE)), aes(label=as.character(snp)),
                           nudge_y  = 9 -max(dat2$PVALUE),
                           box.padding   = 0.5, 
                           point.padding = 0.5,
                           segment.color = 'grey50')
lz <- lz + geom_point(aes(x=dat2$POS[dat2$PVALUE == max(dat2$PVALUE)], y=max(dat2$PVALUE)), 
                      colour="purple", shape=18, size=3) +
  theme(legend.position = c(0.96,0.88))

plot_grid(lz, gene_panel, align="v", axis="l", nrow=2, rel_heights = c(2,1))


dev.off()