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

png(paste0("Manhattan plot of ", title.name,".png"), width=1080, height=540)
manhattan(x=tmp, chr = "scaffold", bp = "base.pair",
		suggestiveline = -log10(0.1/114420), genomewideline = -log10(0.01/114420),
		ylim=c(0,10), cex.main=2,
		snp = "SNP.name", p=var.name, xlab = "Scaffold",
		main=paste("Manhattan plot of", title.name))
dev.off()
}