workDir <- "/home/acari/github/teeth_project"
setwd(workDir)

library(dada2)
library(phyloseq)
library(ggplot2)
library(pheatmap)
library(zCompositions)
library(compositions)
library(selbal)

source("src/mini_compositions_lib.R")

ps <- readRDS("data/ps.rds")

table(tax_table(ps)[, "Phylum"], exclude = NULL)
ps0 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

prevdf <- apply(X = otu_table(ps0),
                MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(ps0),
                     tax_table(ps0))

plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
100-495556*100/3788860
filterPhyla  <- c("Euryarchaeota", "Excavata", "SAR", "Tenericutes", "Saccharibacteria", "Candidate_division_SR1", "Synergistetes")
ps1 <-subset_taxa(ps0, !Phylum %in% filterPhyla)

prevdf1  <- subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps0),color=Phylum)) +
    # Include a guess for parameter
    geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
    scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
    facet_wrap(~Phylum) + theme(legend.position="none")

prevalenceThreshold = 0.05 * nsamples(ps0)
prevalenceThreshold

keepTaxa  <- rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2  <- prune_taxa(keepTaxa, ps1)
PS <- ps2

taxa_names(PS) <- paste0("otu", seq(ntaxa(PS)))

OTU_TABLE <- otu_table(PS)
OTU_TABLE <- as.data.frame(t(OTU_TABLE))
OTU_TABLE <- cbind(OTU = rownames(OTU_TABLE), OTU_TABLE)

TREE <- colnames(otu_table(ps2))

write.csv(OTU_TABLE, "output/OTU_TABLE.csv", row.names = F, quote = F)
write.fasta(as.list(TREE), names = colnames(otu_table(PS)), "output/seqs.fasta")

ps_gen  <- tax_glom(ps2, "Genus", NArm = TRUE)
ps_fam  <- tax_glom(ps2, "Family", NArm = TRUE)

# Heatmap family
taxa_table <- data.frame(seqs = rownames(tax_table(ps_fam)[, "Family"]), family = tax_table(ps_fam)[, "Family"])
rownames(taxa_table) <- 1:nrow(taxa_table)

OTU <- as.data.frame(otu_table(ps_fam))
OTU <- as.data.frame(t(OTU))
OTU <- cbind(seqs = rownames(OTU), OTU)

OTU <- merge(taxa_table, OTU, by = 1)
OTU <- OTU[-1]
OTU <- OTU[-7,]
rownames(OTU) <- OTU$Family
OTU <- OTU[-1]

df.annot <- data.frame(sample = colnames(OTU), group = sample_data(ps_fam)$Group)
rownames(df.annot) <- df.annot$sample
df.annot <- df.annot[-1]
colnames(OTU) <- paste(df.annot$group, colnames(OTU), sep = "")
OTU <- OTU[order(colnames(OTU))]
colnames(OTU) <- sub("C|P|SDP", "", colnames(OTU))

df.annot.col <- list(group = c(C = "red", P="blue", SDP="darkgreen"))
OTU <- OTU[rownames(OTU) != "Mitochondria",]

OTU_2 <- as.data.frame(t(OTU))
OTU_2 <- OTU_2[order(colSums(OTU_2), decreasing = T)]
write.table(OTU_2, "output/OTU.table.txt", sep = "\t", quote = F, row.names = T)

heatmap_family <- pheatmap(log(OTU+10), cluster_cols = F, 
         show_colnames = F, 
         annotation_col = df.annot, annotation_colors = df.annot.col,
         annotation_legend = F,
         gaps_col = c(16,16, 31,31))

svg(filename="figures/heatmap_family.svg", width=10, height=4.5)
heatmap_family
dev.off()

# alpha diversity
ps_alpha_div_shannon <- estimate_richness(ps2, split = TRUE, measure = "Shannon")
ps_alpha_div_simpson <- estimate_richness(ps2, split = TRUE, measure = "Simpson")
ps_alpha_div_chao <- estimate_richness(ps2, split = TRUE, measure = "Chao1")

ps_alpha_div <- cbind(ps_alpha_div_shannon, ps_alpha_div_simpson, ps_alpha_div_chao[1])
ps_alpha_div$group <- sample_data(ps2)$Group
ps_alpha_div$group <- as.character(ps_alpha_div$group)

ps_alpha_div$group[ps_alpha_div$group == "C"] <- "Control"
ps_alpha_div$group[ps_alpha_div$group == "P"] <- "CP"
ps_alpha_div$group[ps_alpha_div$group == "SDP"] <- "CP+T2DM"

shannon_plot <- ggplot(ps_alpha_div, aes(group, Shannon, fill = group))+
    geom_boxplot(alpha = 0.65)+
    ylim(c(2,3.5))+
    theme_classic()+
    ylab("Shannon index")+
    theme(legend.position = "none")+
    scale_fill_manual(values = c("darkred", "darkblue", "darkgreen"))

simpson_plot <- ggplot(ps_alpha_div, aes(group, Simpson, fill = group))+
    geom_boxplot(alpha = 0.65)+
    ylim(c(0.79, 0.99))+
    theme_classic()+
    ylab("Simpson index")+
    theme(legend.position = "none")+
    scale_fill_manual(values = c("darkred", "darkblue", "darkgreen"))

chao1_plot <- ggplot(ps_alpha_div, aes(group, Chao1, fill = group))+
    geom_boxplot(alpha = 0.65)+
    ylim(c(10, 75))+
    theme_classic()+
    ylab("Chao1 index")+
    theme(legend.position = "none")+
    scale_fill_manual(values = c("darkred", "darkblue", "darkgreen"))

svg(filename="figures/shannon_boxplot.svg", width=2.2, height=3, pointsize=12)
shannon_plot
dev.off()

svg(filename="figures/simpson_boxplot.svg", width=2.2, height=3, pointsize=12)
simpson_plot
dev.off()

svg(filename="figures/chao1_boxplot.svg", width=2.2, height=3, pointsize=12)
chao1_plot
dev.off()

A1 <- wilcox.test(ps_alpha_div$Shannon[ps_alpha_div$group == "Control"], ps_alpha_div$Shannon[ps_alpha_div$group == "CP"], alternative = "less")$p.val
B1 <- wilcox.test(ps_alpha_div$Shannon[ps_alpha_div$group == "Control"], ps_alpha_div$Shannon[ps_alpha_div$group == "CP+T2DM"], alternative = "less")$p.val
C1 <- wilcox.test(ps_alpha_div$Shannon[ps_alpha_div$group == "CP"], ps_alpha_div$Shannon[ps_alpha_div$group == "CP+T2DM"], alternative = "less")$p.val

A2 <- wilcox.test(ps_alpha_div$Simpson[ps_alpha_div$group == "Control"], ps_alpha_div$Simpson[ps_alpha_div$group == "CP"], alternative = "less")$p.val
B2 <- wilcox.test(ps_alpha_div$Simpson[ps_alpha_div$group == "Control"], ps_alpha_div$Simpson[ps_alpha_div$group == "CP+T2DM"], alternative = "less")$p.val
C2 <- wilcox.test(ps_alpha_div$Simpson[ps_alpha_div$group == "CP"], ps_alpha_div$Simpson[ps_alpha_div$group == "CP+T2DM"], alternative = "less")$p.val

A3 <- wilcox.test(ps_alpha_div$Chao1[ps_alpha_div$group == "Control"], ps_alpha_div$Chao1[ps_alpha_div$group == "CP"], alternative = "less")$p.val
B3 <- wilcox.test(ps_alpha_div$Chao1[ps_alpha_div$group == "Control"], ps_alpha_div$Chao1[ps_alpha_div$group == "CP+T2DM"], alternative = "less")$p.val
C3 <- wilcox.test(ps_alpha_div$Chao1[ps_alpha_div$group == "CP"], ps_alpha_div$Chao1[ps_alpha_div$group == "CP+T2DM"], alternative = "less")$p.val

wilcox_res <- data.frame(group = c("Control vs CP", "Control vs CP+T2DM", "CP vs CP+T2DM"), shannon = c(A1, B1, C1), simpson = c(A2, B2, C2), chao1 = c(A3, B3, C3))
wilcox_res$shannon <- p.adjust(wilcox_res$shannon, method = "fdr")
wilcox_res$simpson <- p.adjust(wilcox_res$simpson, method = "fdr")
wilcox_res$chao1 <- p.adjust(wilcox_res$chao1, method = "fdr")

write.table(wilcox_res, "output/wilcox_results.txt", quote = F, row.names = F)

# Aitchison MDS
DATA <- as.data.frame(t(OTU))
DATA <- DATA[order(colSums(DATA), decreasing = T)]

Numzz <- colSums(DATA==0)
# Savezz <-  (Numzz<nrow(DATA)*0.48)
Xz <- DATA[1:8]

meta_data <- as.data.frame(sample_data(ps_fam))

TO_SONGBIRD_1 <- Xz[rownames(Xz) %in% rownames(meta_data[meta_data$Group %in% c("SDP", "C")]),]
METADATA_1 <- meta_data[meta_data$Group %in% c("SDP", "C")]

TO_SONGBIRD_2 <- Xz[rownames(Xz) %in% rownames(meta_data[meta_data$Group %in% c("P", "C")]),]
METADATA_2 <- meta_data[meta_data$Group %in% c("P", "C")]

TO_SONGBIRD_3 <- Xz[rownames(Xz) %in% rownames(meta_data[meta_data$Group %in% c("P", "SDP")]),]
METADATA_3 <- meta_data[meta_data$Group %in% c("P", "SDP")]

write.table(t(TO_SONGBIRD_1), "output/TO_SONGBIRD_1.txt", quote = F, sep = "\t")
write.table(METADATA_1, "output/METADATA_1.txt", quote = F, sep = "\t", )

write.table(t(TO_SONGBIRD_2), "output/TO_SONGBIRD_2.txt", quote = F, sep = "\t")
write.table(METADATA_2, "output/METADATA_2.txt", quote = F, sep = "\t", )

write.table(t(TO_SONGBIRD_3), "output/TO_SONGBIRD_3.txt", quote = F, sep = "\t")
write.table(METADATA_3, "output/METADATA_3.txt", quote = F, sep = "\t", )

X <- mcountsprop(X=Xz)

Xclr <-mclr(X)
XAdist <-dist(Xclr)

Xmds <- cmdscale(XAdist ,k=2)
Xmds <- as.data.frame(Xmds)

Xmds <- merge(cbind(rownames(Xmds), Xmds), cbind(rownames(meta_data), meta_data), by = 1)
Xmds$Group <- as.character(Xmds$Group)

Xmds$Group[Xmds$Group == "C"] <- "Control"
Xmds$Group[Xmds$Group == "P"] <- "CP"
Xmds$Group[Xmds$Group == "SDP"] <- "CP+T2DM"

aitchison_mds <- ggplot(Xmds, aes(V1, V2, col = Group))+
    geom_point(size = 2)+
    theme_bw()+
    theme(legend.position = "bottom")+
    scale_color_manual(values = c("red", "blue", "darkgreen"))+
    xlab("MDS1")+
    ylab("MDS2")

svg(filename="figures/aitchison_mds.svg", width=3.5, height=4, pointsize=12)
aitchison_mds
dev.off()

# CoDa dendrogram
V <- mPBclustvar(X)

# ilr coordinates
Xilr <- milr(X,V)

CONTROL <- X[which(!is.na(str_extract(rownames(X), "C"))),]
CP <- X[which(!is.na(str_extract(rownames(X), "P"))),]
CP <- CP[which(is.na(str_extract(rownames(CP), "SDP"))),]
CP_T2DM <- X[which(!is.na(str_extract(rownames(X), "SDP"))),]

svg(filename="figures/coda_dendrogram.svg", width=5, height=3.5)

par(mar = c(10,0,0.5,0))

CoDaDendrogram(acomp(X),V=V,type="lines",
               range=c(-5,5),
               lwd.tree = 1,lwd.leaf = 1, lty.leaf = 2,yaxt="n")
CoDaDendrogram(acomp(CONTROL),V=V,type="lines",
               col="red",lwd.tree = 1,lwd.leaf = 1,yaxt="n",
               lty.leaf = 2,add=TRUE)
CoDaDendrogram(acomp(CP),V=V,type="lines",
               col="blue",lwd.tree = 1,lwd.leaf = 1,yaxt="n",
               lty.leaf = 2,add=TRUE)
CoDaDendrogram(acomp(CP_T2DM),V=V,type="lines",
               col="darkgreen",lwd.tree = 1,lwd.leaf = 1.5,yaxt="n",
               lty.leaf = 2,add=TRUE)

dev.off()