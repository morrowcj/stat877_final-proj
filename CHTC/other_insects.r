rm(list=ls())

args = (commandArgs(trailingOnly=TRUE))

ind <- as.numeric(args[1]) # get the job number to run the corresponding 20 SNPs


load("other_insects.RData")

library(lme4)


cnt <- 1
N <- 20

p.value <- as.data.frame(matrix(NA, nrow=N, ncol=17))
SNP.name <- rep(NA,N)



for (i in 1:20) {

snp.no <- ind*20 + i

for (insect.no in 1:17) {
    lookup <- data.frame(FID=rownames(SNP.add),
                     snp=SNP.add[,snp.no])

    tmp.data2 <- merge(lookup, tmp.data, by = 'FID')
    tmp.data2 <- cbind(tmp.data2,y=1*(tmp.data2[,2+insect.no]>0))

    m <- glmer(y ~ snp  + Volume + ALA.all +
            SLA.all + BBreakDegDay.all + EFNMean.all +
            CTsum + PGsum  + Age
            # Npct.all + Cpct.all # exclude these two variable due to missing
          + (1|survey.event/FID), family = binomial, 
          data=tmp.data2)

    p.value[cnt,insect.no] <- summary(m)$coef[2,4] # p-value of SNP
	SNP.name[cnt] <- names(SNP.add)[snp.no]
	
    # cat(cnt, insect.no,"\n")
  }
  cnt <- cnt + 1

}

p.value <- cbind(SNP.name, p.value)

names(p.value) <- c("SNP.name", paste0("pval.",names(tmp.data2)[3:19]))


write.csv(result,file=paste0("other_insects_",ind+1,".csv"), row.names = FALSE)