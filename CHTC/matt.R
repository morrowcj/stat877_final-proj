.libPaths(c("/workspace/ting/R/library"))

library(qqman)


load("all_insects_combined.RData")

pseudo.table <- read.csv("pseudo-chrom-key.csv",header = TRUE)

no.insect <- NCOL(all_insects_combined)

for (i in 2:19) {
tmp <- merge(pseudo.table, all_insects_combined[,c(1,i)], by = "SNP.name")

var.name <- names(all_insects_combined)[i]

title.name <- sub("\\.", " ", substr(var.name,6, nchar(var.name)))

png(paste0("Manhattan plot of ", title.name,".png"), width=1080, height=540)
manhattan(x=tmp, chr = "scaffold", bp = "base.pair",
		snp = "SNP.name", p=var.name, xlab = "Scaffold",main=paste("Manhattan plot of", title.name))
dev.off()
}