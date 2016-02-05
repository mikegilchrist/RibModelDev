rm(list=ls())
library(ribModel)

with.phi <- TRUE

if (with.phi) {
  genome <- initializeGenomeObject(file = "../data/singleMixture/simulatedOneMix.fasta", expression.file = "../data/singleMixture/simulatedOneMix_simphi.csv")
} else {
  genome <- initializeGenomeObject(file = "../data/realGenomes/Skluyveri.fasta")
}


sphi_init <- (c(1))
numMixtures <- 1
mixDef <- "allUnique"
geneAssignment <- c(rep(1, genome$getGenomeSize(FALSE)))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)
#parameter <- initializeParameterObject(model="ROC", restart.file="10restartFile.rst")


# initialize MCMC object
samples <- 10000
thining <- 10
adaptiveWidth <- 10
mcmc <- initializeMCMCObject(samples=samples, thining=thining, adaptive.width=adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)


# get model object
model <- initializeModelObject(parameter, "ROC", with.phi)

setRestartSettings(mcmc, "restartFile.rst", adaptiveWidth, TRUE)
#run mcmc on genome with parameter using model
system.time(
  runMCMC(mcmc, genome, model, 8)
)


writeParameterObject(parameter, file="InvGTest_Param.Rdat")
writeMCMCObject(mcmc, file="InvGTest_MCMC.Rdat")
