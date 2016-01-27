library(seqinr)

idx <- sort(sample(1:5530, 2000, replace=F))

phi <- read.table("Scereviciae.phi_estm.csv", sep=",", header=T)
phi2000 <- phi[idx,c(1,2)]
colnames(phi2000) <- c("Gene", "phi")
write.table(phi2000, file = "nse2000.phi.csv", sep=",", row.names=F, col.names=T, quote=F)

genome <- read.fasta(file="Scereviciae.fasta")

names <- getName(genome)

geneidx <- names %in% phi2000[,1]
write.fasta(file = "nse2000.fasta", names = names[geneidx], sequences = genome[geneidx])
