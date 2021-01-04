argsv <- commandArgs(trailingOnly = TRUE)

run <- as.numeric(argsv[1]) #run
total <- as.numeric(argsv[2]) #total number of annotations

base <- read.table("./BASE/BF.BASE")
annot <- read.table("SN_data.txt", header=FALSE, stringsAsFactors = FALSE)
base <- base[1,1]


bf <- read.table(paste("SN/BF.",run,sep=""), header=F)
  
bf_annot=bf[1,1]
  
stat= -2*(base - bf_annot)
pval <- pchisq(stat, df=1, lower.tail = F)

if(pval < 0.05/total){
  cat("STOP\n", file="runfile.txt")
}else{
  cat("RUN\n", file="runfile.txt")
}
