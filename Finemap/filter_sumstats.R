argsv <- commandArgs(trailingOnly = TRUE)
chr <- argsv[1]

#procedure:
#subset sumstats and write out
#write out summary file
#read in posterior prob files
#estimate 95% credible set in each and subset again
#write-out updated summary stats

#for base model
#files <- list.files("results")
#files <- files[grep("results", files)]
#files <- files[grep("bin", files)]

#for model with annotations
files <- list.files("FINAL")
files <- files[grep(paste0("chr", chr), files)]
files <- files[grep("Results", files)]

for(file in files){
  pp <- read.table(paste0("FINAL/",file), header=TRUE, stringsAsFactors = FALSE, colClasses = c("character"))
  pp$Posterior_Prob <- as.numeric(pp$Posterior_Prob)
  total95 <- sum(pp$Posterior_Prob)*0.95
  total90 <- sum(pp$Posterior_Prob)*0.90
  pp <- pp[order(pp$Posterior_Prob, decreasing = TRUE),]
  pp$cred95 <- NA
  j <- 0
  for(i in 1:nrow(pp)){
    j <- pp$Posterior_Prob[i] + j
    pp$cred95[i] <- ifelse(j < total95, TRUE, FALSE)
  }
  
  if(file == files[1]){
    all <- pp
  }else{
    all <- rbind(all,pp)
  }
}

perDROP <- round(nrow(all[all$cred95 != TRUE,])/nrow(all)*100)
print(paste0("removing ", perDROP,"%", " of variants"))

final <- all[all$cred95 == TRUE,]


write.table(final, paste0("chr",chr,".pruned_sumstats.txt"), sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE)
