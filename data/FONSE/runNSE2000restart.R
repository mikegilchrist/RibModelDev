library(ribModel)
rm(list=ls())
#read genome

genome.file <- "nse2000.fasta"
phi.file <- "nse2000.phi.csv"
mut.file <- "Scereviciae.mut.csv"
sel.file <- "Scereviciae.sel.csv"
res.file <- "150restartFile.rst"

from.good.values <- FALSE
genome <- initializeGenomeObject(genome.file)

sphi_init <- 1
numMixtures <- 1
mixDef <- "allUnique"
#geneAssignment <- c(rep(1,448), rep(1,513), rep(2,457), rep(1, 3903))
#geneAssignment <- c(rep(1,448), rep(1,513), rep(2,457))
#geneAssignment <- c(rep(1,448), rep(2,457))
#geneAssignment <- c(rep(1,500), rep(2,500))
geneAssignment <- rep(1,length(genome))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, model="FONSE", split.serine = TRUE,
                                       mixture.definition = mixDef, restart.file = res.file)

# initialize MCMC object
samples <- 400 
thinning <- 10 
adaptiveWidth <- 50000 
mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)
# get model object
model <- initializeModelObject(parameter, "FONSE")

if(from.good.values) {
  parameter$initMutationCategories(c(mut.file), 1)
	parameter$initSelectionCategories(c(sel.file), 1)
	phi <- read.table(phi.file, header = T, sep = ",")
	phi.values <- phi[,2]
				  
	parameter$initializeSynthesisRateByList(phi.values)
}

setRestartSettings(mcmc, "restartFile.rst", 100, TRUE)
#run mcmc on genome with parameter using model
system.time(
  runMCMC(mcmc, genome, model, 16)
)

#plots log likelihood trace, possibly other mcmc diagnostics in the future
pdf("NSE2000restartMCMC.pdf")
plot(mcmc)
loglik.trace <- mcmc$getLogLikelihoodTrace()
acf(loglik.trace)

# plots different aspects of trace
trace <- parameter$getTraceObject()
plot(trace, what = "MixtureProbability")
plot(trace, what = "SPhi")
plot(trace, what = "ExpectedPhi")

mixtureAssignment <- unlist(lapply(1:genome$getGenomeSize(),  function(geneIndex){parameter$getEstimatedMixtureAssignmentForGene(samples, geneIndex)}))
expressionValues <- unlist(lapply(1:genome$getGenomeSize(), function(geneIndex){
  expressionCategory <- parameter$getSynthesisRateCategoryForMixture(mixtureAssignment[geneIndex])
  parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(samples, geneIndex, expressionCategory)
}))
expressionValues <- log10(expressionValues)
obs.phi <- log10(read.table(phi.file, sep=",", header=T)[, 2])
cat(length(expressionValues), "\n")
cat(length(obs.phi), "\n")
cat(range(expressionValues, na.rm = T))
plot(NULL, NULL, ylim=range(expressionValues, na.rm = T) + c(-0.1, 0.1), xlim=range(obs.phi) + c(-0.1, 0.1), 
     main = "Synthesis Rate", xlab = "true values", ylab = "estimated values")
upper.panel.plot(obs.phi[mixtureAssignment == 1], expressionValues[mixtureAssignment == 1], col="black")
legend("topleft", legend = paste("Mixture Element", 1:numMixtures), 
       col = ribModel:::.mixtureColors[1:numMixtures], lty = rep(1, numMixtures), bty = "n")

#plot(parameter, what = "Mutation")
#plot(parameter, what = "Selection")
dev.off()


#plot(trace, what = "Expression", geneIndex = 999, mixture = 2)

pdf("NSE2000restartCSP.pdf", width = 11, height = 12)
mixture <- 1
plot(trace, what = "Mutation", mixture = mixture)
plot(trace, what = "Selection", mixture = mixture)
#dev.off()
# plots model fit (cub plot)
#plot(model, genome, parameter, samples = samples*0.1, mixture = mixture, main = "Codon Usage Plot")

names.aa <- aminoAcids()
selection <- c()
mutation <- c()
#csp <- read.table("../ribModel/data/simulated_CSP0.csv", sep=",", header=T)
mut <- read.table(mut.file, sep=",", header=FALSE)
sel <- read.table(sel.file, sep=",", header=FALSE)
#idx.eta <- grepl(pattern = "[A-Z].[A-Z]{3}.Delta.eta", x = as.character(csp[,1]))
#idx.mu <- grepl(pattern = "[A-Z].[A-Z]{3}.log.mu", x = as.character(csp[,1]))
for(aa in names.aa)
{
  if(aa == "M" || aa == "W" || aa == "X") next
  codons <- AAToCodon(aa, T)
  for(i in 1:length(codons))
  {
    selection <- c(selection, parameter$getSelectionPosteriorMeanForCodon(mixture, samples*0.1, codons[i]))
    mutation <- c(mutation, parameter$getMutationPosteriorMeanForCodon(mixture, samples*0.1, codons[i]))
  }
}
#plot(NULL, NULL, xlim=range(csp[idx.mu, 2], na.rm = T), ylim=range(mutation), 
#     main = "Mutation", xlab = "true values", ylab = "estimated values")
plot(NULL, NULL, xlim=range(mut[,3], na.rm = T), ylim=range(mutation), 
     main = "Mutation", xlab = "true values", ylab = "estimated values")
upper.panel.plot(mut[,3], mutation)
#plot(NULL, NULL, xlim=range(csp[idx.eta, 2], na.rm = T), ylim=range(selection), 
#     main = "Selection", xlab = "true values", ylab = "estimated values")
plot(NULL, NULL, xlim=range(sel[,3], na.rm = T), ylim=range(selection), 
     main = "Selection", xlab = "true values", ylab = "estimated values")
upper.panel.plot(sel[,3], selection)
dev.off()
