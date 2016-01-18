rm(list=ls())
library(ribModel)


parameter <- loadParameterObject(c("RFPObject1.Rdat", "RFPObject2.Rdat"))
mcmc <- loadMCMCObject(c("MCMCObject1.Rdat", "MCMCObject2.Rdat"))
genome <- initializeGenomeObject(file = "../data/rfp/rfp.counts.by.codon.and.gene.GSE63789.wt.csv", FALSE)


trace <- parameter$getTraceObject()
samples <- mcmc$getSamples()
pdf("test1.pdf")


plot(mcmc) #plots the whole logliklihood trace

#Here I take a subset of the trace values for the logliklihood trace and plot them.
#The primary reason for doing this is the "jump" that throws the scale of the graph
#at the beginning is removed by taking out the beginning values.
loglik.trace <- mcmc$getLogLikelihoodTrace()
start <- length(loglik.trace) * 0.5 #the multiplier determines how much of the beginning trace is 
#eliminated.

logL <- logL <- mean(loglik.trace[start:length(loglik.trace)]) #get the mean for the subset
plot(loglik.trace[start:length(loglik.trace)], type="l", main=paste("logL:", logL), xlab="Sample", ylab="log(Likelihood)")
grid (NULL,NULL, lty = 6, col = "cornsilk2")


plot(trace, what = "MixtureProbability") #right now, will be straight line (mix =1)
plot(trace, what = "Mphi")
plot(trace, what = "Sphi")
plot(trace, what = "ExpectedPhi")
loglik.trace <- mcmc$getLogLikelihoodTrace()
acf(loglik.trace)
acf(loglik.trace[start:length(loglik.trace)])
dev.off()



pdf("test2.pdf", width = 11, height = 20)
#plot(trace, what = "Expression", geneIndex = 905) #used to make sure gene stabalized, not really needed now
plot(trace, what = "Alpha", mixture = 1)
plot(trace, what = "LambdaPrime", mixture = 1)
dev.off()







pdf("test3.pdf")
cat <- 1
proposal <- FALSE
alphaList <- numeric (61)
lambdaPrimeList <- numeric (61)
waitingTimes <- numeric(61)
phiList <- numeric(genome$getGenomeSize())
ids <- numeric(genome$getGenomeSize())
codonList <- codons()
i <- 1
for (i in 1:61)
{
  codon <- codonList[i]
  alphaList[i] <- parameter$getCodonSpecificPosteriorMean(cat, samples * 0.5, codon, 0)
  lambdaPrimeList[i] <- parameter$getCodonSpecificPosteriorMean(cat, samples * 0.5, codon, 1)
  waitingTimes[i] <- alphaList[i] * lambdaPrimeList[i]
}

for (geneIndex in 1:genome$getGenomeSize()) {
  phiList[geneIndex] <- parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(samples * 0.5, geneIndex, 1)
}

for (i in 1:genome$getGenomeSize())
{
  g <- genome$getGeneByIndex(i, FALSE)
  ids[i] <- g$id
}






plot(NULL, NULL, xlim=range(alphaList, na.rm = T), ylim=range(lambdaPrimeList), 
     main = "Correlation Between Alpha and Lambda Prime", xlab = "alpha", ylab = "lambdaPrime")
upper.panel.plot(alphaList, lambdaPrimeList)

dev.off()
