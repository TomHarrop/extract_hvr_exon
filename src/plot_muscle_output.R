#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(ggtree)
library(ggplot2)

# dist_dt <- fread("output/align_consensus/consensus.dist",
#                  header = FALSE)
# setnames(dist_dt, names(dist_dt), c("V1", dist_dt[, V1]))

# read the MUSCLE phylip tree
nw <- treeio::read.newick(snakemake@input[["tree"]])
nw_rooted <- ape::root(nw, outgroup = "ref")
tree_data <- treeio::as.treedata(nw_rooted)

# turn it into a tree
# d <- as.dist(as.matrix(data.frame(dist_dt, row.names = "V1")))
# nj <- ape::bionj(d)
# nj_rooted <- ape::root(nj, outgroup = "ref")
# 
# tree_data <- treeio::as.treedata(nj_rooted)

# dist dt like ref?
weird_drones <- paste0("BB", c(12, 18, 24, 35, 36, 49, 54))

# weird_dt <- copy(dist_dt)
# setorder(weird_dt, ref)
# weird_dt[, weird := V1 %in% weird_drones]
# weird_dt[, data.table(t(summary(ref))), by = weird]

# cluster and cut into alleles
my_dist <- as.dist(cophenetic(nw))
hc <- hclust(my_dist)
memb <- cutree(hc, k = snakemake@params[["alleles"]]) # arbitrary

# annotate
annot <- data.table(tip_label = nw_rooted$tip.label)
annot[, indiv := gsub("_.*", "", tip_label)]
dup_indivs <- annot[duplicated(annot, by = "indiv"), unique(indiv)]
annot[indiv %in% dup_indivs, tech_dup := TRUE]
annot[!indiv %in% dup_indivs, tech_dup := NA]
annot[indiv %in% weird_drones, deletion := "Low coverage region"]
annot[!indiv %in% weird_drones, deletion := "OK"]
annot[indiv == "ref", deletion := "Reference"]
annot[, cluster := memb[tip_label]]

# plot
gp <- ggtree(tree_data, layout = "circular") %<+% annot +
    theme(legend.position="right",
          text = element_text(size = 8)) +
    scale_colour_viridis_d(guide = guide_legend(title = "Allele")) +
    scale_fill_brewer(palette = "Set1",
                      guide = guide_legend(title = NULL)) +
    geom_tiplab2(mapping = aes(colour = as.factor(cluster)),
                 size = 3,
                 offset = 0.0025) +
    geom_tippoint(mapping = aes(fill = deletion),
                  colour = "white",
                  shape = 21,
                  size = 2)

# write
ggsave(snakemake@output[["plot"]],
       gp,
       device = cairo_pdf,
       width = 10,
       height = 7.5,
       units = "in")

# log
sessionInfo()
