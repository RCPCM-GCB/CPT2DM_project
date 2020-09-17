workDir <- "/home/acari/github/teeth_project"
setwd(workDir)

library(phangorn)
library(seqinr)
library(ggtree)

genomes_taxonomy <- read.csv("data/genomes_taxonomy.csv", stringsAsFactors = F)
genomes_taxonomy <- data.frame(id = genomes_taxonomy$SEQ.ID, taxa = apply(genomes_taxonomy[1:3], 1, function(x) paste(x, collapse = " ")), stringsAsFactors = F)

tree <- read.tree("data/tree.nwk")

mags_table <- read.table("output/mags_table.txt", header = T, stringsAsFactors = F)
true_names <- c(mags_table$mag.name, mags_table$min.bac.name)
true_names <- true_names[true_names != "SEQF1468"]

tree.p <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% true_names])

genomes_taxonomy.sbs <- genomes_taxonomy[as.character(genomes_taxonomy$id) %in% tree.p$tip.label,]

for (k in as.character(genomes_taxonomy.sbs$id)){
    tree.p$tip.label[tree.p$tip.label == k] <- genomes_taxonomy.sbs$taxa[genomes_taxonomy.sbs$id == k]
}

tree_plot <- ggtree(tree.p)+
    geom_tiplab(align=TRUE)+
    xlim(0, 3.5)

svg(filename="figures/tree.svg", width=10, height=10)
tree_plot
dev.off()