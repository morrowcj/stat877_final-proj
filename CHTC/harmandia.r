rm(list=ls())

args = (commandArgs(trailingOnly=TRUE))

ind <- as.numeric(args[1]) # get the job number to run the corresponding 20 SNPs


load("stat877.RData")

library(lme4)

result <- as.data.frame(matrix(NA,nrow=20,ncol=10))

for (i in 1:20) {

snp.no <- ind*20 + i

lookup <- data.frame(FID=rownames(SNP.add),
                     snp=SNP.add[,snp.no])

tmp.data2 <- merge(lookup, tmp.data, by = 'FID')

m <- glmer(y ~ snp  + Volume + ALA.all +
        SLA.all + BBreakDegDay.all + EFNMean.all +
        CTsum + PGsum  + Age
      + (1|FID), family = binomial, data=tmp.data2)

result[i,1]  <- names(SNP.add)[snp.no]
result[i,-1] <-  summary(m)$coef[-1,4]

}			

names(result) <- c("SNP.name", "pval.SNP", "pval.Volume", "pval.ALA.all", 
        "pval.SLA.all", "pval.BBreakDegDay.all", "pval.EFNMean.all",
        "pval.CTsum", "pval.PGsum","pval.Age")

write.csv(result,file=paste0("harmandia_",ind+1,".csv"), row.names = FALSE)