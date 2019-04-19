#!/usr/bin/Rscript

## Title ----
## R script to create manhatten plots and to ID common/significant SNPs

## Arguments/flags ----
# Set default arguments
p.cutoff = .01 # -c ; default = .01
save = FALSE # -s ; default = FALSE
make.plots = TRUE # --no-plot ; default = TRUE
  hist.plot = FALSE # -h ; default = FALSE
  gg.plot = FALSE # -g ; default = FALSE
    interactive.plot = FALSE # --interactive-plots ; default = FALSE
  scaff.plots = FALSE # --scaff-plot ; default = FALSE
  circ.plot = FALSE # --circ ; default = FALSE
print.results = TRUE  # --no-print ; default = TRUE
# change arguments if given via command line
{args = commandArgs()
  if ("-c" %in% args) {
    p.cutoff <- as.numeric(args[which(args == "-c")+1])
    stopifnot(is.numeric(p.cutoff))
    message("p.cutoff = ",p.cutoff)
  }
  if ("-s" %in% args) {save = TRUE} # save objects made here
  if ("--no-plot" %in% args) {make.plots = FALSE} # don't make or save any plots
  if ("-h" %in% args) {hist.plot = TRUE} # make and save histograms
  if ("-g" %in% args) {gg.plot = TRUE} # make the gg manhattan plots
  if ("--interactive-plots" %in% args) {interactive.plot = TRUE} # make int. plots
  if ("--scaff-plot" %in% args) {scaff.plots = TRUE} # make select scaffold plots
  if ("--circ" %in% args) {circ.plot = TRUE} # make circular manhattan plots
  if ("--no-print" %in% args) {print.results = FALSE} # don't print results
}

## Setup ----
# load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(qqman)
  library(CMplot)
  library(qvalue)

})

# load data
data("all_insects_combined")
pseudo.table <- read.csv("data/pseudo-chrom-key.csv",header = TRUE)
# count scaffolds/pseudo-chromosomes
scaffold.n <- length(unique(pseudo.table$scaffold))
# count insects
no.insect <- NCOL(all_insects_combined)
# create empty lists for snps
insect.sig.snps <- list()
running.sig.snps <- NULL
running.sig.scaf <- NULL
# calculate q value matrix (outside loop)
q.matrix <- mapply(all_insects_combined[,-1],
                   FUN = function(x){
                     qvalue(x, fdr.level = .05)$qvalues
                     })
q.df <- data.frame("SNP.name" = as.character(all_insects_combined[,1]),q.matrix)
names(q.df) <- gsub(names(q.df),pattern = "pval",replacement = "qval")
# create things before the loop
gplot.list <- list() # full manhattan plots
filtered.gplot.list <- list() # select scaffold manhattan plots
gg.objects <- NULL

# loop over insects
for (i in 2:no.insect) {
  ## Calculations ----
  # create temparary data
  tmp <- merge(pseudo.table,
               all_insects_combined[,c(1,i)],
               by = "SNP.name")
  # insect variable name (with "pval." prefix)
  var.name <- names(all_insects_combined)[i]
  # add q value and fdr corrections (within loop)
  q.list <- qvalue(tmp[,4],fdr.level = .05)
  tmp$q.vals <- q.list$qvalues # Storey's q
  tmp$lfdr <- q.list$lfdr # local fdr
  tmp$signf <- q.list$significant # q < .05?
  # make a list of significant snps
  significant.snps <- tmp[tmp$signf,]
  ## add relevant columns to the full list
  insect.sig.snps[[var.name]] <- list(
    "SNP" = significant.snps[,"SNP.name"],
    "scaffold" = significant.snps[,"scaffold"],
    "df" = significant.snps) # full data frame just in case
  ## add to the running list of significant SNPs
  running.sig.snps <- unique(c(running.sig.snps,
                              as.character(significant.snps[,"SNP.name"])))
  ## ditto for scaffolds
  running.sig.scaf <- unique(c(running.sig.scaf,
                               as.character(significant.snps[,"scaffold"])))
  # create a subset to plot (p < p.cutoff)
  pval.col <- grep(names(tmp),pattern = "pval", value = TRUE)
  q.cutoff <- max(tmp[tmp[,pval.col] <= p.cutoff,"q.vals"])
  plot.subset <- subset(tmp, tmp[,pval.col] < p.cutoff)
  # Remove "." from insect names to make pot titles and save names
  title.name <- gsub("\\.", " ", substr(var.name,6, nchar(var.name)))
  save.name <- gsub(title.name, pattern = "\\s{1,}",replacement = "-")

  if (print.results | make.plots) {
    # print loop iteration and insect name
    cat(i-1, title.name,"\n")
  }

  ## Plots ----
  if (make.plots) {
    if (hist.plot) {
      # Histograms
      ## If there were no difference in genotype among present/absent groups,
      ## we would expect a histogram of p-values that is relatively flat (uniform).
      ## if there is a poisson-like distribution of p-values (long tail to right)
      ## then there is likely a difference in 'expression' i.e. genotype.
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
    if (gg.plot) {
      # Full manhattan ggplots
      ## load significant SNP objects
      gg.env <- new.env()
      load("data/significant-SNP-list.RData",envir = gg.env)
      ## make list of snps to highlight
      snps.of.interest <- unique(unlist(gg.env$incommon.snps$SNP.name))
      ## prepare data to plot:
      dat.gg <- plot.subset %>%
        # Compute chromosome size
        group_by(scaffold) %>% summarize("scaff.len" = max(base.pair)) %>%
        # Calculate cumulative position of each chromosome
        mutate(tot = cumsum(scaff.len)-scaff.len) %>% select(-scaff.len) %>%
        # Add this info to the initial dataset
        left_join(plot.subset, ., by = c("scaffold" = "scaffold")) %>%
        # Add a cumulative position of each SNP
        arrange(scaffold, base.pair) %>% mutate( BPcum = base.pair + tot) %>%
        # add highlight and annotation info (here same as significant)
        mutate(is_highlight = ifelse(SNP.name %in% snps.of.interest,
                                     yes = TRUE, no = FALSE))
      ### prepare x axis
      axisdf = dat.gg %>% group_by(scaffold) %>%
        summarize(center = (max(BPcum) + min(BPcum))/2)
      ## ggplot of q values
      gplot.list[[save.name]] <-
        ggplot(dat.gg, aes(x = BPcum, y = -log10(q.vals))) +
        # significance lines
        geom_hline(yintercept = -log10(c(.05,.01)),
                   color = c("darkblue","blue"), linetype = "dashed") +
        # snp pointss
        geom_point(aes(color = as.factor(scaffold)), alpha = .8, size = 1.3) +
        scale_color_manual(values = rep(c("darkgrey","black"),scaffold.n)) +
        # custom X and Y axes
        scale_x_continuous(labels = axisdf$scaffold, breaks = axisdf$center) +
        scale_y_continuous(expand = c(0,0), limits = c(-log10(q.cutoff), 3.5)) +
        # add highlighted points and make them slightly bigger:
        geom_point(data = subset(dat.gg, is_highlight),
                   color = "red", size = 1.5) +
        # label plot:
        labs(y = expression({-log[10](symbol(q))}),
             x = "Genome Location (scaffold)",
             title = paste("Manhattan Plot of",title.name)) +
        # custom theme
        theme_bw() +
        theme(legend.position = "none",
              panel.border = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.text.x = element_blank(),
              text = element_text(face = "bold",size = 12))
      # save the plots as images:
      ggsave(filename = paste0("figures/q-manhattan_ggplot/ggManhattan-plot_", save.name,".png"),
             plot = gplot.list[[save.name]],width = 14.4, height = 7.2, units = "in")
      # add to save list on last iteration
      if(i == no.insect){
        gg.objects <- c(gg.objects, "gplot.list")
        }

        if(interactive.plot){
          ## load plotly
          suppressPackageStartupMessages({library(plotly)})
          ### prepare SNP description
          dat.gg$text <- paste("SNP:",dat.gg$SNP.name,
                               "\nPosition:",dat.gg$base.pair,
                               "nScaffold:",dat.gg$scaffold,
                               "\nq-LOD score:",-log10(dat.gg$q.vals)
          )
          ggplotly(gplot.list[[save.name]], tooltip = "text")
        }
    }
    if (scaff.plots) {
      # ggplot of q values, only with certain scaffolds
      ## load sig snp data
      small.env <- new.env()
      load("data/significant-SNP-list.RData",envir = small.env)
      snps.of.interest <- unique(unlist(small.env$incommon.snps$SNP.name))
      filtered.scaffolds <- pseudo.table[match(small.env$incommon.snps$SNP.name,
                                               pseudo.table$SNP.name),]
      ## Prepare plot data
      filt.dat <- tmp %>%
        # remove scaffolds that don't have sig snps for >1 insect
        filter(scaffold %in% filtered.scaffolds$scaffold) %>%
        mutate(scaffold = factor(scaffold))
      filt.dat.gg <- filt.dat %>%
        # Compute chromosome size
        group_by(scaffold) %>% summarize("scaff.len" = max(base.pair)) %>%
        # Calculate cumulative position otf each chromosome
        mutate(tot = cumsum(scaff.len)-scaff.len) %>% select(-scaff.len) %>%
        # Add this info to the initial dataset
        left_join(filt.dat, ., by = c("scaffold" = "scaffold")) %>%
        # Add a cumulative position of each SNP
        arrange(scaffold, base.pair) %>% mutate( BPcum = base.pair + tot) %>%
        # # add highlight and annotation info (here same as significant)
        mutate(is_highlight = ifelse(SNP.name %in% snps.of.interest,
                                     yes = TRUE, no = FALSE))
      ### prepare x axis
      axisdf = filt.dat.gg %>% group_by(scaffold) %>%
        summarize(center = (max(BPcum) + min(BPcum))/2)
      ## create vertical grid for highlighted snps
      sig.lines <- filt.dat.gg %>%
        filter(SNP.name %in% snps.of.interest) %>% select(BPcum)
      ### create plot layer for the lines
      sig.line.layer <- geom_vline(xintercept = unlist(sig.lines),
                                   col = "darkgreen", linetype = "dashed",
                                   size = .05)
      # Make the plot
      filtered.gplot.list[[save.name]] <-
        ggplot(filt.dat.gg, aes(x = BPcum, y = -log10(q.vals))) +
        geom_hline(yintercept = -log10(c(.05,.01)),
                   color = c("darkblue","blue"), size = .5, linetype = "dashed") +
        sig.line.layer + # add vertical line layer
        geom_point(aes(color = as.factor(scaffold)), alpha = .8, size = 1) +
        scale_color_manual(values = rep(c("darkgrey","black"),scaffold.n)) +
        # custom X axis
        scale_x_continuous(labels = axisdf$scaffold, breaks = axisdf$center) +
        scale_y_continuous(expand = c(0,0), limits = c(-log10(1), 3.5)) +
        # add highlighted points:
        geom_point(data = subset(filt.dat.gg, is_highlight),
                   color = "red", size = 1) +
        labs(y = expression({-log[10](symbol(q))}),
             x = "Genome Location (scaffolds with shared signifance only)",
             title = paste("Manhattan Plot of",title.name),
             subtitle = "(select scaffolds)") +
        # custom theme
        theme_bw() +
        theme(legend.position = "none",
              panel.border = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.major.y = element_blank(),
              axis.text.x = element_blank(),
              text = element_text(face = "bold",size = 12))
      # save plots to images
      ggsave(filename = paste0("figures/sig-scaffold-q-manhattan_ggplot/ggManhattan-plot-filtered_",
                               save.name,".png"),
             plot = filtered.gplot.list[[save.name]],
             width = 14.4, height = 7.2, units = "in")
      # add the plot objects to the list on the last iteration
      if(i == no.insect){
        gg.objects <- c(gg.objects, "filtered.gplot.list", "sig.line.layer")
      }

    }

  }
}
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

if (scaff.plots & make.plots) {
  library(gridExtra)
  library(grid)
# alter the filtered plots slightly, to fit them all onto one plot:
  harm <- filtered.gplot.list[["Harmandia"]] +
    labs(title = NULL, subtitle = NULL,
         y = "Harmandia",
         x = NULL)
  phyl <- filtered.gplot.list[["Phyllocolpa"]] +
    labs(title = NULL, subtitle = NULL,
         y = "Phyllocolpa",
         x = NULL)
  pet <- filtered.gplot.list[["Petiole-Gall"]] +
    labs(title = NULL, subtitle = NULL,
         y = "Petiole Gall",
         x = NULL)
  grn.aph <- filtered.gplot.list[["Green-Aphids"]] +
    labs(title = NULL, subtitle = NULL,
         y = "Green Aphid",
         x = NULL)
  smk.aph <- filtered.gplot.list[["Smokey-Aphids"]] +
    labs(title = NULL, subtitle = NULL,
         y = "Smokey Aphid",
         x = NULL)
  # grid.plots <- list(harm,phyl,pet,grn.aph,smk.aph,code =  '# Stack the plots on top of each other
  # multi.plot.scaff <- grid.arrange(ncol = 1,
  #   # plots to stack
  #   harm, phyl, pet, grn.aph, smk.aph,
  #   # shared y axis
  #   left = textGrob(expression({-log[10](symbol(q))}), rot = 90, vjust = .5,
  #                   gp = gpar(fontface = "bold", cex = 1.3)),
  #   # shared x axis
  #   bottom = textGrob("Genome Location (scaffolds with shared signifance only)", hjust = .5,
  #                     gp = gpar(fontface = "bold", cex = 1.2)),
  #   # shared title
  #   top = textGrob("Manhattan Plot of Shared Significant SNPs", hjust = .5,
  #                  gp = gpar(fontface = "bold", cex = 1.5))
  #   )')
  # Stack the plots on top of each other
  multi.plot.scaff <- grid.arrange(ncol = 1,
    # plots to stack
    harm, phyl, pet, grn.aph, smk.aph,
    # shared y axis
    left = textGrob(expression({-log[10](symbol(q))}), rot = 90, vjust = .5,
                    gp = gpar(fontface = "bold", cex = 1.3)),
    # shared x axis
    bottom = textGrob("Genome Location (scaffolds with shared signifance only)", hjust = .5,
                      gp = gpar(fontface = "bold", cex = 1.2)),
    # shared title
    top = textGrob("Manhattan Plot of Shared Significant SNPs", hjust = .5,
                   gp = gpar(fontface = "bold", cex = 1.5))
    )
  ggsave(filename = "figures/gg-multiManhattanPlot.png",
         plot = multi.plot.scaff,
         # width = 14.4, height = 7.2, units = "in"
         width = 16, height = 8, units = "in")

  gg.objects <- c(gg.objects, "multi.plot.scaff", "grid.plots")

}

## Results ----
# matrix of significant
sig.mat <- matrix(NA, ncol = length(insect.sig.snps), nrow = length(running.sig.snps))
colnames(sig.mat) <- names(insect.sig.snps)
rownames(sig.mat) <- running.sig.snps


for (j in colnames(sig.mat)) {
  sig.mat[,j] <- ifelse(running.sig.snps %in% insect.sig.snps[[j]][["SNP"]], yes = 1, no = 0)
}
## count number of insect associations per SNP
sig.df <- as.data.frame(sig.mat)
sig.df$n.assoc <- rowSums(sig.mat)
sig.df$SNP.name <-  rownames(sig.mat)
# sig.df <- merge(sig.df,q.df, by = "SNP.name") # add in q values
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
# total sig.snps by insect
sig.n.per.insect <- lapply(insect.sig.snps,
                           function(x){
                             return(c("SNPs" = length(unique(x[[1]])),
                                      "scaffolds" = length(unique(x[[2]]))))
                           })

sig.n.table <- data.frame(NULL)
for (i in names(sig.n.per.insect)) {
  table.row <- gsub(gsub(i, pattern = "pval\\.", replacement = ""),
                    pattern = "\\.", replacement = " ")
  sig.n.table[table.row,"SNPs"] <-  sig.n.per.insect[[i]][1]
  sig.n.table[table.row,"scaffolds"] <- sig.n.per.insect[[i]][2]
}
sig.n.table <- sig.n.table[order(sig.n.table$SNPs, decreasing = TRUE),]

scaffolds.uniq <- unique(filtered.scaffolds$scaffold.name)
assoc.combos.uniq <- unique(snp.insects)
if (print.results) {
  cat("Total Unique Significant SNPs:\n")
  print(length(running.sig.snps))

  cat("Total Unique Significant scaffolds:\n")
  print(length(running.sig.scaf))

  cat("Total Significant SNPS and Scaffolds per Insect:\n")
  print(sig.n.per.insect)

  cat("Significant Associations:\n")
  print(snp.insects)

  cat("Unique Association Pairs:\n")
  print(assoc.combos.uniq)

  cat("Unique Significant Scaffolds")
  print(scaffolds.uniq)
}

## Save Results ----
if (save) {
  message("saving results objects ...")
  save(list = c(
    # from "Calculations" section
    "insect.sig.snps","running.sig.snps","running.sig.scaf","snp.insects", "q.df",
    # from "Results" section
    "sig.n.per.insect","sig.df","incommon.snps","filtered.scaffolds",
    "sig.n.table"),
     file = "data/significant-SNP-list.RData")
  # save gg plot objects
  if (!is.null(gg.objects)){
    message("saving ggManhattan plot objects")
    save(list = gg.objects, file = "data/ggManhattan-objects.RData")
  }


}

## End of File ----
message("End")
