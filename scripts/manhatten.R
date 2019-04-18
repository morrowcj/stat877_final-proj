#!/usr/bin/Rscript

##----
## R script to create manhatten plots and to ID common/significant SNPs
##----

##----
# Arguments/flags #
p.cutoff = .01 # -c
save = FALSE # -s
make.plots = TRUE # --no-plot
  man.plot = FALSE # -m
  hist.plot = FALSE # -h
  gg.plot = FALSE # -g
    interactive.plot = FALSE # --interactive-plots
  small.plots = FALSE # --scaff-plot
  circ.plot = FALSE # --circ
print.results = TRUE  # --no-print

args = commandArgs()
if ("-c" %in% args) {
  p.cutoff <- as.numeric(args[which(args == "-c")+1])
  stopifnot(is.numeric(p.cutoff))
  message("p.cutoff = ",p.cutoff)
}
if ("-s" %in% args) {save = TRUE}
if ("--no-plot" %in% args) {make.plots = FALSE}
if ("-m" %in% args) {man.plot = TRUE}
if ("-h" %in% args) {hist.plot = TRUE}
if ("-g" %in% args) {gg.plot = TRUE}
if ("--interactive-plots" %in% args) {interactive.plot = TRUE}
if ("--scaff-plot" %in% args) {small.plots = TRUE}
if ("--circ" %in% args) {circ.plot = TRUE}
if ("--no-print" %in% args) {print.results = FALSE}
# if ("-h" %in% args | "--help" %in% args) {cat("")}
##----

##----
# Setup #
## load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(qqman)
  library(CMplot)
  library(qvalue)

})

## load data
data("all_insects_combined")
pseudo.table <- read.csv("data/pseudo-chrom-key.csv",header = TRUE)

## count scaffolds/pseudo-chromosomes
scaffold.n <- length(unique(pseudo.table$scaffold))
## count insects
no.insect <- NCOL(all_insects_combined)

## create empty lists for snps
insect.sig.snps <- list()
running.sig.snps <- NULL
running.sig.scaf <- NULL


## calculate q value matrix
q.matrix <- mapply(all_insects_combined[,-1], FUN = function(x){
  qvalue(x, fdr.level = .05)$qvalues
})
q.df <- data.frame("SNP.name" = as.character(all_insects_combined[,1]),q.matrix)
names(q.df) <- gsub(names(q.df),pattern = "pval",replacement = "qval")


## create things before the loop
gplot.list <- list()
filtered.gplot.list <- list()
##----

##----
# Body
## loop over insects
for (i in 2:no.insect) {
  ## create temparary data
  tmp <- merge(pseudo.table,
               all_insects_combined[,c(1,i)],
               by = "SNP.name")
  ## insect name
  var.name <- names(all_insects_combined)[i]
  ## add q value and fdr corrections
  q.list <- qvalue(tmp[,4],fdr.level = .05)
  tmp$q.vals <- q.list$qvalues # Storey's q
  tmp$lfdr <- q.list$lfdr # local fdr
  tmp$signf <- q.list$significant # q < .05?
  significant.snps <- tmp[tmp$signf,]
  insect.sig.snps[[var.name]] <- list("SNP" = significant.snps[,"SNP.name"],
                                      "scaffold" = significant.snps[,"scaffold"],
                                      "df" = significant.snps)
  running.sig.snps <- unique(c(running.sig.snps,
                              as.character(significant.snps[,"SNP.name"])))
  running.sig.scaf <- unique(c(running.sig.scaf,
                               as.character(significant.snps[,"scaffold"])))

  # create a subset to plot (p < p.cutoff)
  pval.col <- grep(names(tmp),pattern = "pval", value = TRUE)
  q.cutoff <- max(tmp[tmp[,pval.col] <= p.cutoff,"q.vals"])
  plot.subset <- subset(tmp, tmp[,pval.col] < p.cutoff)

  # Remove . for title
  title.name <- gsub("\\.", " ", substr(var.name,6, nchar(var.name)))
  save.name <- gsub(title.name, pattern = "\\s{1,}",replacement = "-")

  if (print.results | make.plots) {
    cat(i-1, title.name,"\n")
  }


  ##----
  # Plots #
  if (make.plots) {
    # Histograms
      ## If there were no difference in genotype among present/absent groups,
      ## we would expect a histogram of p-values that is relatively flat (uniform).
      ## if there is a poisson-like distribution of p-values (long tail to right)
      ## then there is likely a difference in 'expression' i.e. genotype.
    if (hist.plot) {
      # p value:
      png(paste0("figures/p-hist/p-histogram-of_", save.name, ".png"), width=720, height=720)
      hist(tmp[,var.name], probability = TRUE,
           xlab = "p value",
           main = paste("Histogram of", title.name, "p values"),
           cex.main = 2)
      dev.off()
      # q value:
      png(paste0("figures/q-hist/q-histogram-of_", save.name, ".png"), width=720, height=720)
      hist(tmp[,"q.vals"], probability = TRUE,
           xlab = "q value",
           main = paste("Histogram of", title.name, "q values"),
           cex.main = 2)
      dev.off()
    }
    # Manhattan Plots
    if (man.plot) { # qqman version
      # bonf. plot: # After Bonferroni correction, add two lines for 0.1 and 0.01 FWER
      png(paste0("figures/p-manhattan/Manhattan-plot_", save.name,".png"),
          width = 1080, height = 540)
      manhattan(x = plot.subset, chr = "scaffold", bp = "base.pair",
      # bonf correction
      suggestiveline = -log10(0.1/114420), genomewideline = -log10(0.01/114420),
      ylim = c(-log10(p.cutoff),10), cex.main = 2,
      snp = "SNP.name", p = var.name, xlab = "Scaffold",
      main = paste("Manhattan plot of", title.name))
      dev.off()
      # q value
  		png(paste0("figures/q-manhattan/Manhattan-plot_", save.name,".png"),
  		    width = 1080, height = 540)
  		manhattan(x = plot.subset, chr = "scaffold", bp = "base.pair",
      		suggestiveline = -log10(0.05), genomewideline = -log10(0.01),
      		ylim = c(-log10(q.cutoff),3.5), cex.main = 2,
      		snp = "SNP.name", p = "q.vals", xlab = "Scaffold",
      		main = paste("Manhattan plot of", title.name),
      		ylab = expression({-log[10](symbol(q))}))
      dev.off()
    }
    if (gg.plot) { # ggplot version
      ## load plotly
      suppressPackageStartupMessages({library(plotly)})
      ## load interesting SNP list
      gg.env <- new.env()
      load("data/significant-SNP-list.RData",envir = gg.env)
      snps.of.interest <- unique(unlist(gg.env$incommon.snps))
      ## prepare data:
      dat.gg <- plot.subset %>%
        # filter(SNP.name %in% snps.of.interest) %>%
        # Compute chromosome size
        group_by(scaffold) %>%
        summarize("scaff.len" = max(base.pair)) %>%
        # Calculate cumulative position of each chromosome
        mutate(tot = cumsum(scaff.len)-scaff.len) %>%
        select(-scaff.len) %>%
        # Add this info to the initial dataset
        left_join(plot.subset, ., by = c("scaffold" = "scaffold")) %>%
        # Add a cumulative position of each SNP
        arrange(scaffold, base.pair) %>%
        mutate( BPcum = base.pair + tot) %>%
        # # add highlight and annotation info (here same as significant)
        mutate(is_highlight = ifelse(SNP.name %in% snps.of.interest,
                                     yes = TRUE, no = FALSE))

      ### prepare x axis
      axisdf = dat.gg %>% group_by(scaffold) %>%
        summarize(center = (max(BPcum) + min(BPcum))/2)

      ### prepare SNP description
      dat.gg$text <- paste("SNP:",dat.gg$SNP.name,
                           "\nPosition:",dat.gg$base.pair,
                           "nScaffold:",dat.gg$scaffold,
                           "\nq-LOD score:",-log10(dat.gg$q.vals)
                           )

      ## ggplot of q values
      gplot.list[[save.name]] <- ggplot(dat.gg, aes(x = BPcum, y = -log10(q.vals))) +
        geom_hline(yintercept = -log10(c(.05,.01)), color = c("blue","red"), linetype = "dashed") +
        geom_point(aes(color = as.factor(scaffold)), alpha = .8, size = 1.3) +
        scale_color_manual(values = rep(c("darkgrey","black"),scaffold.n)) +
        # custom X axis
        scale_x_continuous(labels = axisdf$scaffold, breaks = axisdf$center) +
        scale_y_continuous(expand = c(0,0)) +
        ylim(-log10(q.cutoff), 3.5) +
        # add highlighted points:
        geom_point(data = subset(dat.gg, is_highlight), color = "orange", size = 1.2) +
        labs(y = expression({-log[10](symbol(q))}),
             x = "Genome Location (scaffold)",
             title = paste("Manhattan Plot of",title.name)) +
        # custom theme
        #theme_bw()
        theme(legend.position = "none",
              panel.border = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.text.x = element_blank(),
              text = element_text(face = "bold",size = 12))

      ggsave(filename = paste0("figures/q-manhattan_ggplot/ggManhattan-plot_", save.name,".png"),
             plot = gplot.list[[save.name]],width = 14.4, height = 7.2, units = "in")
        if(!interactive.plot){
          ggplotly(gplot.list[[save.name]], tooltip = "text")
        }
    }
    if (small.plots) {
      ## ggplot of q values, only with certain scaffolds
      small.env <- new.env()
      load("data/significant-SNP-list.RData",envir = small.env)
      snps.of.interest <- unique(unlist(small.env$incommon.snps))
      filtered.scaffolds <- pseudo.table[match(small.env$incommon.snps$SNP.name,pseudo.table$SNP.name),]
      # filt.dat.gg <- dat.gg %>% filter(scaffold %in% filtered.scaffolds$scaffold)

      filt.dat <- tmp %>% filter(scaffold %in% filtered.scaffolds$scaffold) %>%
        mutate(scaffold = factor(scaffold))
      filt.dat.gg <- filt.dat %>%
        # Compute chromosome size
        group_by(scaffold) %>%
        summarize("scaff.len" = max(base.pair)) %>%
        # Calculate cumulative position otf each chromosome
        mutate(tot = cumsum(scaff.len)-scaff.len) %>%
        select(-scaff.len) %>%
        # Add this info to the initial dataset
        left_join(filt.dat, ., by = c("scaffold" = "scaffold")) %>%
        # Add a cumulative position of each SNP
        arrange(scaffold, base.pair) %>%
        mutate( BPcum = base.pair + tot) %>%
        # # add highlight and annotation info (here same as significant)
        mutate(is_highlight = ifelse(SNP.name %in% snps.of.interest,
                                     yes = TRUE, no = FALSE))
      # filt.dat.gg$scaffold <- factor(filt.dat.gg$scaffold, levels = unique(filt.dat.gg$scaffold))
      ### prepare x axis
      axisdf = filt.dat.gg %>% group_by(scaffold) %>%
        summarize(center = (max(BPcum) + min(BPcum))/2)


      filtered.gplot.list[[save.name]] <- ggplot(filt.dat.gg, aes(x = BPcum, y = -log10(q.vals))) +
        geom_hline(yintercept = -log10(c(.05,.01)), color = c("blue","red"), linetype = "dashed") +
        geom_point(aes(color = as.factor(scaffold)), alpha = .8, size = 1.3) +
        scale_color_manual(values = rep(c("darkgrey","black"),scaffold.n)) +
        # custom X axis
        scale_x_continuous(labels = axisdf$scaffold, breaks = axisdf$center) +
        scale_y_continuous(expand = c(0,0)) +
        ylim(-log10(1), 3.5) +
        # add highlighted points:
        geom_point(data = subset(filt.dat.gg, is_highlight), color = "orange", size = 1.2) +
        labs(y = expression({-log[10](symbol(q))}),
             x = "Genome Location (scaffold)",
             title = paste("Manhattan Plot of",title.name),subtitle = "(select scaffolds)") +
        # custom theme
        #theme_bw()
        theme(legend.position = "none",
              panel.border = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.text.x = element_blank(),
              text = element_text(face = "bold",size = 12))
      # filtered.gplot.list[[save.name]]
      ggsave(filename = paste0("figures/sig-scaffold-q-manhattan_ggplot/ggManhattan-plot-filtered_",
                               save.name,".png"),
             plot = filtered.gplot.list[[save.name]],
             width = 14.4, height = 7.2, units = "in")
    }

  }
  ##----
}
##----
# multi-tracked manhattan plots
if(circ.plot & make.plots){
  ## Prepare data
  plot.dat <- merge(pseudo.table,
                    all_insects_combined,
                    by = "SNP.name")
  names(plot.dat) <- gsub(names(plot.dat),pattern = "pval\\.",replacement = "")
  filtered.plot.dat <- plot.dat %>% select(SNP.name,scaffold,base.pair,
                                           Harmandia,Phyllocolpa,Petiole.Gall,
                                           Green.Aphids,Smokey.Aphids)

  suppressPackageStartupMessages({library("CMplot")})
  setwd("figures/CMplots")
  # few insects
  CMplot(Pmap = filtered.plot.dat, multracks = TRUE,
         col = c("grey30","grey60"), LOG10 = TRUE,
         threshold = c(.05/114420,.01/114420),threshold.lty = c(1,2),
         threshold.col = c("black","grey"),
         amplify = TRUE, signal.col = c("red","green"), signal.cex = c(1,1),
         ylim = c(-log10(.0001),10),
         file = "jpg", dpi = 300, file.output = TRUE, verbose = TRUE
  )


  setwd("../..")
  stopifnot(getwd() == "/home/morrowcj/Documents/school/stats/877/final-project")
}
##----

##----
# Results #
if (print.results) {
# total sig.snps by insect
lapply(insect.sig.snps,
       function(x){
         snps.n = length(x$SNP)
         scaf.n = length(x$scaffold)
         return(c("SNPs" = snps.n, "scaffolds" = scaf.n))
        })

  # total sig.snps from all insects
  length(running.sig.snps)
}


sig.mat <- matrix(NA, ncol = length(insect.sig.snps), nrow = length(running.sig.snps))
colnames(sig.mat) <- names(insect.sig.snps)
rownames(sig.mat) <- running.sig.snps


for (j in colnames(sig.mat)) {
  sig.mat[,j] <- ifelse(running.sig.snps %in% insect.sig.snps[[j]][["SNP"]], yes = 1, no = 0)
}
## count number of insect associations per SNP
sig.df <- as.data.frame(cbind(sig.mat,"n.assoc" = rowSums(sig.mat)))
sig.df$SNP.name <-  rownames(sig.mat)
incommon.snps <- sig.df %>% filter(n.assoc > 1) #84

snp.insects <- apply(X = incommon.snps %>%
        select(names(insect.sig.snps)),
      MARGIN = 1,
      FUN = function(x){
        gsub(names(x)[x == 1],
             pattern = "pval\\.",
             replacement = "")
        })
names(snp.insects) <- incommon.snps$SNP.name

filtered.scaffolds <- pseudo.table[match(incommon.snps$SNP.name,pseudo.table$SNP.name),]
filtered.scaffolds$scaffold.name <- gsub(filtered.scaffolds$SNP.name,
                                         pattern =  ":\\d*",replacement = "")

if (print.results) {
  cat("Total Unique Significant SNPs:\n")
  print(length(running.sig.snps))

  cat("Total Unique Significant scaffolds:\n")
  print(length(running.sig.scaf))

  cat("Total Significant SNPS and Scaffolds per Insect:\n")
  print(sig.n.per.insect <- lapply(insect.sig.snps,
               function(x){
                 return(c("SNPs" = length(unique(x[[1]])),
                          "scaffolds" = length(unique(x[[2]]))))
               }))

  cat("Significant Associations:\n")
  print(snp.insects)

  cat("Unique Association Pairs:\n")
  snp.insects %>% unique() %>% print()

  cat("Unique Significant Scaffolds")
  unique(filtered.scaffolds$scaffold.name)

}

if (save) {
  save(list = c("insect.sig.snps","running.sig.snps","snp.insects",
                "sig.n.per.insect","sig.df","incommon.snps","filtered.scaffolds","filtered.gplot.list","gplot.list"),
     file = "data/significant-SNP-list.RData")
}
##----

# example of what qvalue()$significant means:
  # manhattan(x = significant.snps, chr = "scaffold", bp = "base.pair",
  #           suggestiveline = -log10(0.05),
  #           genomewideline = -log10(0.01),
  #           ylim = c(-log10(q.cutoff),3), cex.main = 2,
  #           snp = "SNP.name", p="q.vals", xlab = "Scaffold",
  #           main=paste("Manhattan plot of", title.name))

message("End")
