.libPaths(c("/workspace/ting/R/library"))

library(qqman)


load("all_insects_combined.RData")

pseudo.table <- read.csv("pseudo-chrom-key.csv",header = TRUE)

no.insect <- NCOL(all_insects_combined)

for (i in 2:19) {
tmp <- merge(pseudo.table, all_insects_combined[,c(1,i)], by = "SNP.name")

var.name <- names(all_insects_combined)[i]

# Remove . for title 
title.name <- gsub("\\.", " ", substr(var.name,6, nchar(var.name)))

cat(i, title.name,"\n")

# After Bonferroni correction, add two lines for 0.1 and 0.01 FWER

png(paste0("QQ plot of ", title.name,".png"), width=720, height=720)
	
qq(tmp[,4], cex.main=2, main=paste("QQ plot of", title.name))
		
dev.off()
}