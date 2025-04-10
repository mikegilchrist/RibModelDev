library(ribModel)
rm(list=ls())
#read genome
genome <- initializeGenomeObject(file = "../data/twoMixtures/simulatedAllUniqueR.fasta")

#initialize parameter object
sphi_init <- 2
numMixtures <- 2
mixDef <- "allUnique"
geneAssignment <- c(rep(1,500), rep(2,500))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)


#init from true values
phiVals <- parameter$readPhiValues( "../data/twoMixtures/simulatedSelectionSharedR_phi.csv" )
parameter$initializeSynthesisRateByRandom(phiVals)
parameter$initMutationCategories(c("../data/twoMixtures/simulated_mutation0.csv", "../data/twoMixtures/simulated_mutation1.csv"), 2)
parameter$initSelectionCategories(c("../data/twoMixtures/simulated_selection0.csv", "../data/twoMixtures/simulated_selection1.csv") , 2)


# initialize MCMC object
samples <- 100
thinning <- 10
adaptiveWidth <- 10
mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth, 
                     est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)
# get model object
model <- initializeModelObject(parameter, "ROC")

setRestartSettings(mcmc, "restartFile.rst", adaptiveWidth, TRUE)
#run mcmc on genome with parameter using model
system.time(
  runMCMC(mcmc, genome, model)
)
writeParameterToCSV(parameter,"MutationMixture1.csv", CSP="Mutation", 1)
writeParameterToCSV(parameter,filename="MutationMixture2.csv", CSP="Mutation", 2)
writeParameterToCSV(parameter,filename="SelectionMixture1.csv", CSP="Selection", 1)
writeParameterToCSV(parameter,filename="SelectionMixture2.csv", CSP="Selection", 2)
#plots log likelihood trace, possibly other mcmc diagnostics in the future
plot(mcmc)

# plots different aspects of trace
trace <- parameter$getTraceObject()
pdf("simulated_Genome_allUnique_startCSP_True_startPhi_true_adaptSphi_True.pdf")
plot(trace, what = "MixtureProbability")
plot(trace, what = "SPhi")
plot(trace, what = "ExpectedPhi")
dev.off()
plot(trace, what = "Expression", geneIndex = 905)
pdf("simulated_Genome_allUnique_startCSP_True_startPhi_true_adaptSphi_True_mix1.pdf", width = 11, height = 12)
plot(trace, what = "Mutation", mixture = 1)
plot(trace, what = "Selection", mixture = 1)

# plots model fit (cub plot)
plot(model, genome, samples = samples*0.1, mixture = 1, main = "S. kluyveri Chr (A,B,Cleft) Codon Usage Plot")
dev.off()

plot(parameter, what = "Mutation", main = "Mutation Correlation, Not shared")
plot(parameter, what = "Selection", main = "Selecion Correlation, Not shared")

##-----------------------------------------##

#creates a plot Mixture Prob vs log10(phi) (Not part of the package)
num.genes <- length(genome)
mixtureAssignment <- do.call("rbind", lapply(1:num.genes,  function(geneIndex){parameter$getEstimatedMixtureAssignmentProbabilitiesForGene(samples*0.1, geneIndex)}))
plot(NULL, NULL, xlim = c(0,1), ylim = c(-2, 1), xlab = "Probability beeing in Mixture 2", ylab = "log10(phi)")
colors <- c("black", "red")
for(mixture in 1:2)
{
  genes.in.mixture <- which(round(mixtureAssignment[, mixture]+1) == mixture)
  expressionCategory <- parameter$getSynthesisRateCategoryForMixture(mixture)
  
  # need expression values to know range
  num.genes <- length(genes.in.mixture)
  expressionValues <- unlist(lapply(genes.in.mixture, function(geneIndex){
    parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(samples*0.1, geneIndex, expressionCategory)
  }))
  points(mixtureAssignment[genes.in.mixture, mixture], log10(expressionValues), col=colors[mixture])
}
