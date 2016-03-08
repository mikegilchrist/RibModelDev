rm(list=ls())
library(ribModel)

parameter <- loadParameterObject(c("ROCParameter.Rdat"))
mcmc <- loadMCMCObject(c("MCMCObject.Rdat"))