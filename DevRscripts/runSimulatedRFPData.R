rm(list=ls())
library(ribModel)



#read genome
genome <- initializeGenomeObject(file = "../data/rfp/simulatedRFPData1.csv", FALSE)


#initialize parameter object
sphi_init <- c(2)
numMixtures <- 1
mixDef <- "allUnique"
geneAssignment <- c(rep(1, genome$getGenomeSize(F)))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, model= "RFP", split.serine = TRUE, mixture.definition = mixDef)
#parameter <- new(RFPParameter, "10000restartFile.rst")
parameter <- initializeParameterObject(model = 'RFP', restart.file = "10000restartFile.rst")

#init from "true" values

phiVals <- parameter$readPhiValues( "../data/rfp/RFPPhiValues.csv" )
parameter$initializeSynthesisRateByList(phiVals)
parameter$initMutationSelectionCategories(c("../data/rfp/RFPAlphaValues.csv"), 1, "Alpha")
parameter$initMutationSelectionCategories(c("../data/rfp/RFPLambdaPrimeValues.csv"), 1, "LambdaPrime")


# initialize MCMC object
samples <- 100
thining <- 10
adaptiveWidth <- 10
mcmc <- initializeMCMCObject(samples=samples, thining=thining, adaptive.width=adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)
# get model object
model <- initializeModelObject(parameter, "RFP")

setRestartSettings(mcmc, "restartFile.rst", adaptiveWidth, TRUE)

#run mcmc on genome with parameter using model
system.time(
  runMCMC(mcmc, genome, model, 8)
)


# plots different aspects of trace
trace <- parameter$getTraceObject()


pdf("RFP_Genome_allUnique_startCSP_True_startPhi_true_adaptSphi_True.pdf")


plot(mcmc) #plots the whole logliklihood trace

#Here I take a subset of the trace values for the logliklihood trace and plot them.
#The primary reason for doing this is the "jump" that throws the scale of the graph
#at the beginning is removed by taking out the beginning values.
loglik.trace <- mcmc$getLogLikelihoodTrace()
start <- length(loglik.trace) * 0.5 #the multiplier determines how much of the beginning trace is 
#eliminated.

logL <- logL <- mean(loglik.trace[start:length(loglik.trace)]) #get the mean for the subset
plot(tail(loglik.trace, length(loglik.trace)/2), type="l", main=paste("logL:", logL), xlab="Sample", ylab="log(Likelihood)", xaxt="n")
axis(1, xaxp=c(floor(length(loglik.trace)/2), length(loglik.trace), 5), las=2)
grid (NULL,NULL, lty = 6, col = "cornsilk2")


plot(trace, what = "MixtureProbability") #right now, will be straight line (mix =1)
plot(trace, what = "Mphi")
plot(trace, what = "Sphi")
plot(trace, what = "ExpectedPhi")
loglik.trace <- mcmc$getLogLikelihoodTrace()
acf(loglik.trace)
acf(loglik.trace[start:length(loglik.trace)])
dev.off()



pdf("RFP_CSP_Values_Mixture1.pdf", width = 11, height = 20)
#plot(trace, what = "Expression", geneIndex = 905) #used to make sure gene stabalized, not really needed now
plot(trace, what = "Alpha", mixture = 1)
plot(trace, what = "LambdaPrime", mixture = 1)
dev.off()







pdf("ConfidenceIntervalsForAlphaAndLambdaPrime.pdf")
cat <- 1
proposal <- FALSE
alphaList <- numeric (61)
lambdaPrimeList <- numeric (61)
waitingTimes <- numeric(61)
alpha.ci <- matrix(0, ncol=2, nrow=61)
lambdaPrime.ci <- matrix(0, ncol=2, nrow=61)
phiList <- numeric(genome$getGenomeSize(F))
ids <- numeric(genome$getGenomeSize(F))
codonList <- codons()
for (i in 1:61)
{
  codon <- codonList[i]
  alphaList[i] <- parameter$getCodonSpecificPosteriorMean(cat, samples * 0.5, codon, 0)
  alphaTrace <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1, codon, 0)
  alpha.ci[i,] <- quantile(alphaTrace[(samples * 0.5):samples], probs = c(0.025,0.975))
  
  
  lambdaPrimeList[i] <- parameter$getCodonSpecificPosteriorMean(cat, samples * 0.5, codon, 1)
  lambdaPrimeTrace <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1, codon, 1)
  lambdaPrime.ci[i,] <- quantile(lambdaPrimeTrace[(samples * 0.5):samples], probs = c(0.025,0.975))
  waitingTimes[i] <- alphaList[i] * lambdaPrimeList[i]
}
plot(NULL, NULL, xlim=range(1:61, na.rm = T), ylim=range(alpha.ci), 
     main = "Confidence Intervals for Alpha Parameter", xlab = "Codons", 
     ylab = "Estimated values", axes=F) 
confidenceInterval.plot(x = 1:61, y = alphaList, sd.y = alpha.ci)
axis(2)
axis(1, tck = 0.02, labels = codonList[1:61], at=1:61, las=2, cex.axis=.6)

plot(NULL, NULL, xlim=range(1:61, na.rm = T), ylim=range(lambdaPrime.ci), 
     main = "Confidence Intervals for LambdaPrime Parameter", xlab = "Codons", 
     ylab = "Estimated values", axes=F) 
confidenceInterval.plot(x = 1:61, y = lambdaPrimeList, sd.y = lambdaPrime.ci)
axis(2)
axis(1, tck = 0.02, labels = codonList[1:61], at=1:61, las=2, cex.axis=.6)

for (geneIndex in 1:genome$getGenomeSize(F)) {
  phiList[geneIndex] <- parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(samples * 0.5, geneIndex, 1)
}

for (i in 1:genome$getGenomeSize(F))
{
  g <- genome$getGeneByIndex(i, FALSE)
  ids[i] <- g$id
}




#corrolation between RFPModel and Premal's data
X <- read.table("../data/rfp/codon.specific.translation.rates.table.csv", header = TRUE, sep =",")
X <- X[order(X[,1]) , ]

XM <- matrix(c(X[,1], X[,2]), ncol = 2, byrow = FALSE)
Y <- data.frame(codonList[-c(62,63,64)], waitingTimes)
colnames(Y) <- c("Codon", "PausingTime")
Y <- Y[order(Y[,1]) , ]

plot(NULL, NULL, xlim=range(XM[,2], na.rm = T), ylim=range(Y[,2]), 
     main = "Correlation Between Premal and RFP Model Pausing Times", xlab = "True Values", ylab = "Run Values")
upper.panel.plot(XM[,2], Y[,2])
dev.off()



trueAlphaValues <- read.csv("../data/rfp/RFPAlphaValues.csv")
trueLambdaPrimeValues <- read.csv("../data/rfp/RFPLambdaPrimeValues.csv")
truePhiValues <- read.csv( "../data/rfp/RFPPhiValues.csv")

trueAlphaValues <- trueAlphaValues[order(trueAlphaValues[,1]) , ]
trueLambdaPrimeValues <- trueLambdaPrimeValues[order(trueLambdaPrimeValues[,1]) , ]


lambdaPrimeDF <- data.frame(codonList[-c(62,63,64)], lambdaPrimeList)
colnames(lambdaPrimeDF) <- c("Codon", "LambdaPrimeValue")
lambdaPrimeDF <- lambdaPrimeDF[order(lambdaPrimeDF[,1]) , ]
plot(trueLambdaPrimeValues[,2], lambdaPrimeDF[,2], xlim=range(trueLambdaPrimeValues[,2]),
     ylim=range(lambdaPrimeDF[,2]))
upper.panel.plot(trueLambdaPrimeValues[,2], lambdaPrimeDF[,2])

