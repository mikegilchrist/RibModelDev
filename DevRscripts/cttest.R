rm(list=ls())
library(ribModel)
#read genome

with.phi <- FALSE

genome <- initializeGenomeObject(file = "../data/realGenomes/Skluyveri.fasta")

sphi <- (c(1,1,1))
numMixtures <- 3
mixDef <- "allUnique"
geneAssignment <- c(rep(1,448), rep(1,513), rep(2,457), rep(1, 3403), rep(3, 500))

## defaults from initializeParameterObject
expressionValues = NULL
model = "ROC"
split.serine = TRUE
mixture.definition = "allUnique"
mixture.definition.matrix = NULL
restart.file = NULL
mutation_prior_sd = 0.35

if(is.null(restart.file)){
  if(length(sphi) != numMixtures){
    stop("Not all mixtures have an Sphi value assigned!\n")
  }
  
  if(length(genome) != length(geneAssignment)){
    stop("Not all Genes have a mixture assignment!\n")
  }
  
  if(max(geneAssignment) > numMixtures){
    stop("Gene is assigned to non existing mixture!\n")
  }
  
  #TODO: should we check integraty of other values, such as numMixtures being
  #positive?
}

#end initializeParameterObject

if(is.null(mixture.definition.matrix)){ 
  # keyword constructor
  parameter <- new(ROCParameter, as.vector(sphi), numMixtures, geneAssignment, 
                   split.serine, mixture.definition)
}else{
  #matrix constructor
  mixture.definition <- c(mixture.definition.matrix[, 1], 
                          mixture.definition.matrix[, 2])
  parameter <- new(ROCParameter, as.vector(sphi), numMixtures, geneAssignment, 
                   mixture.definition, split.serine)
}


# initialize expression values
if(is.null(expressionValues)){
  parameter$initializeSynthesisRateByGenome(genome, mean(sphi))
}else{
  parameter$initializeSynthesisRateByList(expressionValues)
}

parameter$mutation_prior_sd <- (mutation_prior_sd)

# end initializeROCObject

numMutationCategory <- parameter$numMutationCategories
numSelectionCategory <- parameter$numSelectionCategories

phi <- parameter$getCurrentSynthesisRateForMixture(1) # phi values are all the same initially
#ct <- getInstance()
#names.aa <- ct$getGroupList()
names.aa <- genome$getGroupList()

codonCounts <- getCodonCountsForAA(aa, genome)
numCodons <- dim(codonCounts)[2] - 1



aa <- "A"
ct
qq <- ct$getAAToCodonIndexMapWithoutReference()
ct$AAToAAIndex(aa)
ct$AAToCodon(aa, T)
ct$AAToCodonRange(aa, T)

#-----------------------------------------
# TODO WORKS CURRENTLY ONLY FOR ALLUNIQUE!
#-----------------------------------------
covmat <- vector("list", numMixtures)
mixElement <- 1
for(mixElement in 1:numMixtures){    
  idx <- geneAssignment == mixElement
  csp <- getCSPbyLogit(codonCounts[idx, ], phi[idx])
  
  parameter$initMutation(csp$coef.mat[1,], mixElement, aa)
  parameter$initSelection(csp$coef.mat[2,], mixElement, aa)
}

parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)