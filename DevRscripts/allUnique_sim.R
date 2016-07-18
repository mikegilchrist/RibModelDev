rm(list=ls())
library(ribModel)
seeds <- read.table(file = "seed.txt")[,1]
task.id <- 1 #as.numeric(Sys.getenv("SGE_TASK_ID"))

set.seed(446141)
with.phi <- TRUE

if (with.phi) {
  genome <- initializeGenomeObject(file = "../data/twoMixtures/simulatedAllUniqueR.fasta", 
                                   expression.file = "../data/twoMixtures/simulatedAllUniqueR_phi_withPhiSet.csv") 
} else {
  genome <- initializeGenomeObject(file = "../data/twoMixtures/simulatedAllUniqueR.fasta") 
}

sphi_init <- c(1,1)
numMixtures <- 2
mixDef <- "allUnique"
geneAssignment <- sample(c(1,2), size = length(genome), replace = TRUE, prob = c(0.3, 0.7)) #c(rep(1,500), rep(2,500))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)
#parameter <- initializeParameterObject(restart.file = "5001_simulated_allUnique.rst")

samples <- 10
thining <- 10
adaptiveWidth <- 10
divergence.iteration <- 0
mcmc <- initializeMCMCObject(samples, thining, adaptive.width=adaptiveWidth, est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE) 

setRestartSettings(mcmc, paste(task.id, "_simulated_allUnique.rst", sep=""), adaptiveWidth*50, TRUE) 

model <- initializeModelObject(parameter, "ROC", with.phi = with.phi) 

start <- Sys.time() 
system.time(runMCMC(mcmc, genome, model, 4, divergence.iteration)) 
end <- Sys.time() 
end - start 


pdf(paste(task.id, "_simulated_global_allUnique.pdf", sep=""), width = 11, height = 12) 
plot(mcmc)
loglik.trace <- mcmc$getLogLikelihoodTrace()[-1]
acf(loglik.trace) 
# Currently has an error
#convergence.test(mcmc, n.samples = 500, plot=T) 

trace <- parameter$getTraceObject()
plot(trace, what = "MixtureProbability")
plot(trace, what = "Sphi")
plot(trace, what = "Mphi") 
if (with.phi) { 
  plot(trace, what = "Aphi")
  plot(trace, what = "Sepsilon") 
} 
plot(trace, what = "ExpectedPhi")

mixtureAssignment <- getMixtureAssignmentEstimate(parameter, length(genome), samples*0.1)
expressionValues <- getExpressionEstimatesForMixture(parameter, length(genome), mixtureAssignment, samples*0.1)
expressionValues <- log10(expressionValues) 
obs.phi <- log10(read.table("../data/twoMixtures/simulatedAllUniqueR_phi.csv", sep=",", header=T)[, 2]) 
plot(NULL, NULL, xlim=range(obs.phi) + c(-0.1, 0.1), ylim=range(expressionValues, na.rm = T) + c(-0.1, 0.1), 
     main = "Synthesis Rate", xlab = "True values", ylab = "Estimated values") 
for (k in 1:numMixtures) {
  if (sum(mixtureAssignment == k) == 0) next
  upper.panel.plot(obs.phi[mixtureAssignment == k], expressionValues[mixtureAssignment == k], col=ribModel:::.mixtureColors[k]) 
} 
legend("topleft", legend = paste("Mixture Element", 1:numMixtures), 
       col = ribModel:::.mixtureColors[1:numMixtures], lty = rep(1, numMixtures), bty = "n") 
plot(parameter, what = "Mutation", samples = samples*0.1)
plot(parameter, what = "Selection", samples = samples*0.1) 
dev.off() 

observed.mutation <- c("../data/twoMixtures/simulated_mutation0.csv", "../data/twoMixtures/simulated_mutation1.csv")
observed.selection <- c("../data/twoMixtures/simulated_selection0.csv", "../data/twoMixtures/simulated_selection1.csv")
for (k in 1:numMixtures) {
  if (sum(mixtureAssignment == k) == 0) next
  pdf(paste(task.id, "_simulated_mixture_", k, "_allUnique.pdf", sep=""), width = 11, height = 12) 
  mixture <- k
  plot(trace, what = "Mutation", mixture = mixture)
  plot(trace, what = "Selection", mixture = mixture)
  plot(model, genome, samples = samples*0.1, mixture = mixture, main = "Codon Usage Plot")
  names.aa <- aminoAcids() 
  
  selection.sd <- vector("list")
  mutation.sd <- vector("list")
  
  selection <- c() 
  mutation <- c() 
  codon.storage <- c() 
  csp.m <- read.table(observed.mutation[mixture], sep=",", header=T) 
  csp.e <- read.table(observed.selection[mixture], sep=",", header=T) 
  csp <- rbind(csp.m,csp.e) 
  idx.eta <- 41:80 
  idx.mu <- 1:40 
  for (aa in names.aa) 
  { 
    if (aa == "M" || aa == "W" || aa == "X") next 
    codons <- AAToCodon(aa, T) 
    codon.storage <- c(codon.storage,codons) 
    for (i in 1:length(codons)) 
    {
      # Push back these values to form vectors of length(codons) size.
      mutation <- c(mutation, parameter$getCodonSpecificPosteriorMean(mixture, samples*0.1, codons[i], 0, TRUE))
      mutation.sd <- c(mutation.sd, parameter$getCodonSpecificQuantile(mixture, samples*0.1, codons[i], 0, c(0.025, 0.975), TRUE)) 
      
      selection <- c(selection, parameter$getCodonSpecificPosteriorMean(mixture, samples*0.1, codons[i], 1, TRUE))
      selection.sd <- c(selection.sd, parameter$getCodonSpecificQuantile(mixture, samples*0.1, codons[i], 1, c(0.025, 0.975), TRUE))      
    } 
  }
  mutation.sd <- matrix(unlist(mutation.sd), nrow = 2)
  selection.sd <- matrix(unlist(selection.sd), nrow = 2)

  plot(NULL, NULL, xlim=range(csp[idx.mu, 3], na.rm = T), ylim=range(mutation), 
       main = "Mutation", xlab = "True values", ylab = "Estimated values") 
  upper.panel.plot(x = csp[idx.mu, 3], y = mutation, sd.y = mutation.sd)
  plot(NULL, NULL, xlim=range(csp[idx.eta, 3], na.rm = T), ylim=range(selection), 
       main = "Selection", xlab = "True values", ylab = "Estimated values") 
  upper.panel.plot(x = csp[idx.eta, 3], y = selection, sd.y = selection.sd)
  dev.off() 
} 
