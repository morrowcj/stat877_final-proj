# identify scaffold from SNP name
harm <- read.csv("result/harmandia_pvalue.csv", header = TRUE)

split.list <- strsplit(x = as.character(harm$SNP.name),split = ":")
split.mat <- matrix(unlist(split.list),ncol = 2, byrow = TRUE)
split.mat <- cbind(as.character(harm$SNP.name),split.mat)
colnames(split.mat) <- c("SNP.name","scaffold","base.pair")
split.mat[,"scaffold"] <- gsub(split.mat[,"scaffold"],
                               pattern = "Potra",
                               replacement = "")
split.mat <- split.mat[order(as.numeric(split.mat[,"scaffold"])),]
split.mat[,"scaffold"] <- factor(split.mat[,"scaffold"])
split.df <- as.data.frame(split.mat)

write.csv(split.df,file = "data/pseudo-chrom-key.csv",
          quote = FALSE,row.names = FALSE)



harm.new <- merge(split.df,harm,by = "SNP.name")

# library(ggplot2)
# ggplot(harm.new,aes(x = base.pair, y = pval.SNP, group = scaffold))+
#   geom_bar(stat = "identity") + facet_wrap(~scaffold)
