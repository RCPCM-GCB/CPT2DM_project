workDir <- "/data11/bio/jeniaole/teeth"
setwd("/data11/bio/jeniaole/teeth")

library(phyloseq)
library(phangorn)
library(dada2)
library(ggplot2)
library(DESeq2)
library(DECIPHER)

ref_fasta <- "data/silva_nr_v123_train_set.fa.gz"

miseq_path <- file.path("data", "raw")
filt_path <- file.path("data", "filtered")

fns <- sort(list.files(miseq_path, full.names = TRUE))

fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]

ii <- sample(length(fnFs), 3)
for(i in ii) { print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd")) }
for(i in ii) { print(plotQualityProfile(fnRs[i]) + ggtitle("Rev")) }

filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))

for(i in seq_along(fnFs)) {
    fastqPairedFilter(c(fnFs[[i]], fnRs[[i]]),
                      c(filtFs[[i]], filtRs[[i]]),
                      trimLeft=10, truncLen=c(245, 200),
                      maxN=0, maxEE=2, truncQ=2,
                      compress=TRUE)
}

derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
sam.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepFs) <- sam.names
names(derepRs) <- sam.names

ddF <- dada(derepFs[1:40], err=NULL, selfConsist=T, multithread = T)
ddR <- dada(derepRs[1:40], err=NULL, selfConsist=T, multithread = T)

plotErrors(ddF)
plotErrors(ddR)

dadaFs <- dada(derepFs, err=ddF[[1]]$err_out, pool=T, multithread = T)
dadaRs <- dada(derepRs, err=ddR[[1]]$err_out, pool=T, multithread = T)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

seqtab.all <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
seqtab <- removeBimeraDenovo(seqtab.all)

taxtab <- assignTaxonomy(seqtab, refFasta = ref_fasta)
colnames(taxtab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

seqs <- getSequences(seqtab)
names(seqs) <- seqs 
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)
fit = pml(treeNJ, data=phang.align)

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

group <- gsub(paste0(as.character(0:9), collapse = "|"), "", sam.names)
sample_data <- data.frame(Group = group)
rownames(sample_data) <- sam.names

ps <- phyloseq(tax_table(taxtab),
               otu_table(seqtab, taxa_are_rows = FALSE), 
               phy_tree(fitGTR$tree), sample_data(sample_data))


saveRDS(ps, "data/ps.rds")
saveRDS(seqtab, "data/seqtab.rds")