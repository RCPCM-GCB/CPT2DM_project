workDir <- "/home/acari/github/teeth_project"
setwd(workDir)

library(ggplot2)
library(pheatmap)

CPT2DM_vs_Control <- read.csv("output/differentials_1.tsv", sep = "\t")[-2]
CP_vs_Control <- read.csv("output/differentials_2.tsv", sep = "\t")[-2]
CPT2DM_vs_CP <- read.csv("output/differentials_3.tsv", sep = "\t")[-2]

colnames(CP_vs_Control)[2] <- "effect_size"
colnames(CPT2DM_vs_Control)[2] <- "effect_size"
colnames(CPT2DM_vs_CP)[2] <- "effect_size"

CP_vs_Control$comp_group <- "CP_vs_Control"
CPT2DM_vs_Control$comp_group <- "CPT2DM_vs_Control"
CPT2DM_vs_CP$comp_group <- "CPT2DM_vs_CP"

df <- rbind(CP_vs_Control, CPT2DM_vs_Control, CPT2DM_vs_CP)

df$color[df$effect_size > 0.5] <- "increase"
df$color[df$effect_size < -0.5] <- "dicrease"
df$color[is.na(df$color)] <- "non-significant"

main_barplot <- ggplot(df, aes(reorder(featureid, effect_size), effect_size, fill = color))+
    geom_bar(stat = "identity", width = 0.55)+
    geom_hline(yintercept = 0.5, col = "red", alpha = 0.65)+
    geom_hline(yintercept = -0.5, col = "red", alpha = 0.65)+
    facet_wrap(~comp_group)+
    coord_flip()+
    theme_classic()+
    xlab("Family")+
    ylab("Effect size")+
    scale_fill_manual(values = c("red", "blue", "gray30"))

svg(filename="figures/songbird_barplot.svg", width=10, height=2.5)
main_barplot
dev.off()