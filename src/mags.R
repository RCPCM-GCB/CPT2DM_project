library("ggtree")
library("seqinr")

alignment <- read.alignment("data/alignment.fa", format = "fasta")

D <- dist.alignment(alignment)
D[is.na(D)] <- 1

mags.D <- as.data.frame(as.matrix(D))

mags.sub.D <- mags.D[!is.na(str_extract(colnames(mags.D), "control|CP"))]
mags.sub.D <- mags.sub.D[is.na(str_extract(rownames(mags.sub.D), "control|CP")) ,]

DF <- NULL
for (i in 1:ncol(mags.sub.D)){
    sub.D <- mags.sub.D[i]
    sub.V <- unlist(sub.D)
    names(sub.V) <- rownames(sub.D)
    
    min.bac.name <- names(sub.V[sub.V == min(sub.V)])
    min.bac.dist <- min(sub.V)
    mag.name <- colnames(mags.sub.D)[i]
    
    df.d <- data.frame(mag.name, min.bac.name, min.bac.dist)
    DF <- rbind(DF, df.d)
}

DF.C <- NULL
for (k in as.character(unique(DF$mag.name))){
    DF.C <- rbind(DF.C, DF[DF$mag.name == k,][1,]    )
}

write.table(DF.C, "output/mags_table.txt", sep = "\t", quote = F, row.names = F)
