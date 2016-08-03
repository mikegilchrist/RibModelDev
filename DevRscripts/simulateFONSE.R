library(ribModel)
rm(list=ls())

genome <- initializeGenomeObject(file.path("..", "data", "singleMixture", "genome_2000.fasta"))
mutation.file <- file.path("..", "data", "FONSE", "S.cer.mut.ref.csv")
selection.file <- file.path("..", "data", "FONSE", "selection3ref.csv")
phi.file <- file.path("..", "data", "FONSE", "genome_2000.phi.csv")
genome.out.file <- file.path("..", "data", "FONSE", "fonse3.fasta")

sphi_init <- 1.2
numMixtures <- 1
mixDef <- "allUnique"

geneAssignment <- rep(1, length(genome))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, model= "FONSE", split.serine = TRUE,
                                       mixture.definition = mixDef)
model <- initializeModelObject(parameter, "FONSE")

parameter$initMutationCategories(c(mutation.file), 1)
parameter$initSelectionCategories(c(selection.file), 1)

#phi <- read.table(phi.file, header = T, sep = ",")
#phi.vals <- phi[,2]

phi.vals <- parameter$readPhiValues(phi.file)

parameter$initializeSynthesisRateByList(phi.vals)
model$simulateGenome(genome)

genome$writeFasta(genome.out.file, TRUE)
