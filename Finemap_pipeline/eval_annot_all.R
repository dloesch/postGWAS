#annotation file to check correlations
argsv <- commandArgs(trailingOnly = TRUE)
annot.file <- argsv[1]


base <- read.table("./BASE/BF.BASE")
annot <- read.table("SN_data.txt", header=FALSE, stringsAsFactors = FALSE)
base <- base[1,1]

results <- as.data.frame(1:57)
colnames(results) <- "run"
results$annot <- annot$V1
results$enrich <- NA
results$bf <- NA
results$pvalue <- NA

N <- 57
for (i in 1:N){
  enrich <- read.table(paste("SN/Enrich.",i,sep=""), header=T)
  bf <- read.table(paste("SN/BF.",i,sep=""), header=F)
  results$enrich[i] <- enrich[1,2]
  results$bf[i] <- bf[1,1]
  bf_annot=bf[1,1]
  
  stat= -2*(base - bf_annot)
  results$pvalue[i] <- pchisq(stat, df=1, lower.tail = F)
}


write.table(results, "SN_annot_results.txt", sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE)

results <- results[order(results$pvalue),]
top <- results$annot[1:3] # only keep top 3 annotations

#check correlations
annot <- read.table(annot.file, header=TRUE)

#reformat to match column names
foo <- gsub("-", ".", top)

annot <- annot[colnames(annot) %in% foo]

cor.annot <- as.data.frame(foo)
colnames(cor.annot)[1] <- "annot"
cor.annot$name <- top

for(i in foo){
  for(j in foo)
    cor.annot[[i]][cor.annot$annot == j] <- cor(annot[[i]], annot[[j]])
}

#drop annotations correlated with top annotation
drop <- cor.annot$name[cor.annot[[foo[1]]] > 0.3]
drop <- drop[drop != top[1]]

top <- top[!top %in% drop]
write.table(top, "SN_top_results.txt", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)
