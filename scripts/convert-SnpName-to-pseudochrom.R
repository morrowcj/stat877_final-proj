# identify scaffold from SNP name
## read in harmandia data
harm <- read.csv("result/harmandia_pvalue.csv", header = TRUE)
## split by ":"
split.list <- strsplit(x = as.character(harm$SNP.name),split = ":")
## collapse into a matrix
split.mat <- matrix(unlist(split.list),ncol = 2, byrow = TRUE)
## add in the original SNP name
split.mat <- cbind(as.character(harm$SNP.name),split.mat)
## give the matrix column names
colnames(split.mat) <- c("SNP.name","scaffold","base.pair")
## remove the "Potra"
split.mat[,"scaffold"] <- gsub(split.mat[,"scaffold"],
                               pattern = "Potra",
                               replacement = "")
## Order the matrix according to the scaffold number:
split.mat <- split.mat[order(as.numeric(split.mat[,"scaffold"])),]
## make the scaffold column a factor
split.mat[,"scaffold"] <- factor(split.mat[,"scaffold"])
## turn it into a data frame
split.df <- as.data.frame(split.mat)

# write the key to a file
write.csv(split.df,file = "data/pseudo-chrom-key.csv",
          quote = FALSE,row.names = FALSE)

## To merge in the key
harm <- read.csv("result/harmandia_pvalue.csv", header = TRUE)
split.table <- read.csv("data/pseudo-chrom-key.csv",header = TRUE)
harm.new <- merge(split.table,harm,by = "SNP.name")

library(qqman)
# harm.new$base.pair <- as.numeric(as.character(harm.new$base.pair))
manhattan(head(harm.new,n = 2000),
          chr = "scaffold",
          bp = "base.pair",
          snp = "SNP.name",
          p = "pval.SNP",
          xlab = "Psuedo-chromosome")
