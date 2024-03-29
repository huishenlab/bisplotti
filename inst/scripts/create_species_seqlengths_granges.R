# Create the GRanges objects for each supported species. The goal is to create a lightweight GRanges object for
# seqlengths rather than having to keep larger packages (e.g., BSgenome) as dependencies, speeding up the build system.

#hg19
library(Homo.sapiens)

hg19.gr <- as(seqinfo(Homo.sapiens), "GRanges")

#hg38
library(BSgenome.Hsapiens.UCSC.hg38)

hg38.gr <- as(seqinfo(BSgenome.Hsapiens.UCSC.hg38), "GRanges")

#mm9
library(BSgenome.Mmusculus.UCSC.mm9)

mm9.gr <- as(seqinfo(BSgenome.Mmusculus.UCSC.mm9), "GRanges")

#mm10
library(Mus.musculus)

mm10.gr <- as(seqinfo(Mus.musculus), "GRanges")

#save
save(hg19.gr, file = "bisplotti/data/hg19.gr.rda")
save(hg38.gr, file = "bisplotti/data/hg38.gr.rda")
save(mm9.gr, file = "bisplotti/data/mm9.gr.rda")
save(mm10.gr, file = "bisplotti/data/mm10.gr.rda")
