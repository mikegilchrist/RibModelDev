rm(list=ls())
library(ribModel)


parameter <- loadParameterObject(c("FONSEObject1.Rdat"))
mcmc <- loadMCMCObject(c("MCMCObject1.Rdat", "MCMCObject2.Rdat"))
genome <- initializeGenomeObject(file = "../data/rfp/rfp.counts.by.codon.and.gene.GSE63789.wt.csv", FALSE)
