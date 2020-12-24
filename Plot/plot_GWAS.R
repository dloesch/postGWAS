#!/usr/bin/Rscript

library(data.table)
library(argparser)
library(tools)

#set argument parser
argp <- arg_parser("plot locuszoom-style plot with local LD data")

#reguired arguments
argp <- add_argument(argp, "--results", help="GWAS summary statistics", type = "character")

#optional arguments
argp <- add_argument(argp, "--prefix", help="output prefix", type = "character", default="output")
argp <- add_argument(argp, "--snpcol", help="column name for snp IDs", type = "character", default = NA)
argp <- add_argument(argp, "--pvalcol", help="column name for p-values",type = "character", default = NA)
argp <- add_argument(argp, "--chromcol", help="column name for chromosome", type = "character", default = NA)
argp <- add_argument(argp, "--poscol", help="column name for bp positions", type = "character", default = NA)
argp <- add_argument(argp, "--filetype", help="PDF or PNG for saving plots", type="character", default="PNG")
argp <- add_argument(argp, "--annot", help="add annotations to GWAS-significant peak", type="logical", default=TRUE)
argp <- add_argument(argp, "--annotcol", help="column for adding annotations", type="character", default="snpid")
argp <- add_argument(argp, "--thin", help="thin sub-significant GWAS data", type="logical", default=TRUE)

#parse arguments
argsv <- parse_args(argp)

results.file <- argsv$results

if(is.na(results.file)){
  
  stop("Please provide GWAS summary stats")
}

#read in GWAS data
if(file_ext(results.file) == "RData"){
  load(results.file)
}else{
  data <- fread(results.file, header=TRUE, data.table = FALSE)
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

snpid <- argsv$snpcol
if(is.na(snpid)){
  snpid <- colnames(data)[colnames(data) %in% c("SNP", "snp", "Snp", "variant", "VARIANT", "ID")]
  if(length(snpid) != 1){
    stop("unable to infer SNP ID column name; please provide")
  }
}



data[[chrom]] <- as.numeric(data[[chrom]])

if(argsv$thin == TRUE){
  
  thresh <- 1E-02
  top <- data[data[[pval]] < thresh,]
  temp <- data[data[[pval]] >= thresh,]
  
  chr <- 1
  
  foo <- temp[temp[[chrom]] == chr,]
  
  #keeps 10% of variants below threshold
  thinned <- foo[foo[[snpid]] %in% sample(foo[[snpid]],round(0.1*nrow(foo))),]
  
  for (chr in 2:22){
    foo <- temp[temp[[chrom]] == chr,]
    foo <- foo[foo[[snpid]] %in% sample(foo[[snpid]],round(0.1*nrow(foo))),]
    thinned <- rbind(thinned,foo)
  }
  
  results <- rbind(top, thinned)
  results <- results[order(results[[chrom]], results[[pos]]),]
  results$index <- 1:nrow(results)
  results$logp <- -1*log10(results[[pval]])
  
}else{
 results <- data
 results$index <- 1:nrow(results)
 results$logp <- -1*log10(results[[pval]])
}

gc()

#plot paramgers
y_max <- max(results$logp)
x_max <- max(results$index)
colors <- rep(c("dodgerblue4","skyblue2"),times=11)


#initialize plot
bins <- c()
chr <- 1
i <- 1

if(argsv$filetype == "PNG"){
  png(file=paste0(argsv$prefix, ".mannhattan.png"), width=1000, height=700)
}else{
  pdf(file=paste0(argsv$prefix, ".mannhattan.pdf"),paper = "a4r")
}
#first chromosome
pvals <- results$logp[results[[chrom]] == chr]
index <- results$index[results[[chrom]] == chr]
bins <- c(bins,median(index))
plot(pvals~index, pch=20, col=colors[i], cex=0.8, xaxt="n", ylab="-log10 p-value", 
     xlab="chromosome", main="GWAS results", ylim=c(0, y_max+1), xlim = c(0, x_max+0.5))
i <- i+1

#add remaining chromosomes
for(chr in 2:22){
  pvals <- results$logp[results[[chrom]] == chr]
  index <- results$index[results[[chrom]] == chr]
  bins <- c(bins,median(index))
  points(pvals~index, pch=20, col=colors[i], cex=0.8)
  
  i <- i+1
}

gc()

axis(1, at=bins, labels=c(1:22), las=2)
abline(h= -log10(5E-08),lty = 2)

if(argsv$annot == TRUE){
  for(chr in 1:22){
    pvals <-results$logp[results[[chrom]] == chr]
    if(max(pvals) > -log10(5E-08)){
      annotcol <- argsv$annotcol
      top <- results[[annotcol]][results[[chrom]] == chr & results$logp == max(pvals)][1]
      x <- results$index[results$logp == max(pvals) & results[[chrom]] == chr]
      y <- max(pvals)
      text(x,y+0.25,labels = top, cex = 0.6)
    }
  }
}



dev.off()

#qqplot

#qqplot
#prepare observed pvalues
pvals <- data[[pval]]
obs <- sort(pvals)
obs <- obs[!is.na(obs)]
obs <- obs[is.finite(-log10(obs))]

#prepare expected pvals
exp <- c(1:length(obs))  / (length(1:length(obs))+1) 

x <- exp
y <- obs

if (argsv$thin == TRUE){
  quant.subsample <- function(y, m=100, e=1) {
    # m: size of a systematic sample
    # e: number of extreme values at either end to use
    x <- sort(y)
    n <- length(x)
    quants <- (1 + sin(1:m / (m+1) * pi - pi/2))/2
    sort(c(x[1:e], quantile(x, probs=quants), x[(n+1-e):n]))
    # Returns m + 2*e sorted values from the EDF of y
  }
  n.x <- n.y <- length(obs)
  m <- .001 * max(n.x, n.y)
  e <- floor(0.0005 * max(n.x, n.y))
  x <- quant.subsample(exp, m, e)
  y <- quant.subsample(obs, m, e)
}

if(argsv$filetype == "PNG"){
  png(file=paste0(argsv$prefix, ".qqplot.png"))
}else{
  pdf(file=paste0(argsv$prefix, ".qqplot.pdf"))
}
plot(-log10(x), -log10(y),
     pch=".", cex=4,
     xlab="Expected p-values", ylab="Observed p-values", main="QQ Plot: LARGE-PD GWAS")

abline(0,1,col="red",lwd=3, lty=3)

#add GC lambda to plot
z = qnorm(pvals/2)
lambda = round(median(z^2) / 0.454, 3)
print(lambda)

legend('topleft', legend=paste0("GC Lambda: ",lambda))
dev.off()
