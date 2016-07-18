library(ribModel)
library(MCMCpack)


#TEST randGamma
pdf("randGammaTest.pdf")
shape <- 4
rate <- 2
x <- rgamma(1000, shape = shape, rate = rate)
plot(density(x))
for(i in 1:20)
{
 # abline(v=randGamma(shape, rate))
 # abline(v=rgamma(1, shape=shape, rate=rate), col="green")
#  abline(v=rgamma(1, shape=shape, scale=1/rate), col="red")
  abline(v=sqrt(1/rinvgamma(1, shape=shape, scale=1/rate)), col="purple")
}
dev.off()



shape=2
scale=5
invg <- rinvgamma(n = 1000, shape = shape, scale = scale)
gam <- rgamma(n = 1000, shape = shape, scale = scale)
plot(density(gam))
lines(density(1/invg), col="red")

###############################
#####Inverse Gamma testing#####
###############################

#InvG scale = scale
#G scale = rate, 1/sqrt
#result - seems spot on
shape=6
scale=5
rate=1/scale
invg <- rinvgamma(n = 1000, shape = shape, scale = scale)
plot(density(invg))
for(i in 1:20){
  abline(v=1/sqrt(rgamma(1, shape=shape, scale=rate)))
  abline(v=1/sqrt(randGamma(shape, 1/rate)), col="red")
}


#InvG scale = rate (meaning we could just have a small scale)
#G scale = rate, 1/sqrt
#result - no good
shape=5
scale=4
rate=1/scale
invg <- rinvgamma(n = 1000, shape = shape, scale = rate)
plot(density(invg))
for(i in 1:20){
  abline(v=1/sqrt(rgamma(1, shape=shape, scale=rate)))
}


#InvG scale = rate 
#G scale = scale, 1/sqrt
#result - no good on tail
shape=5
scale=4
rate=1/scale
invg <- rinvgamma(n = 1000, shape = shape, scale = rate)
plot(density(invg))
for(i in 1:20){
  abline(v=1/sqrt(rgamma(1, shape=shape, scale=scale)))
}



#InvG scale = rate 
#G scale = scale, 1/()
#result - also good
shape=6
scale=5
rate=1/scale
invg <- rinvgamma(n = 1000, shape = shape, scale = rate)
plot(density(invg))
for(i in 1:20){
  abline(v=1/(rgamma(1, shape=shape, scale=scale)))
  abline(v=1/(randGamma(shape, rate)), col="red")
}


#InvG scale = scale 
#G scale = scale, 1/()
#result - behind the curve
shape=5
scale=4
rate=1/scale
invg <- rinvgamma(n = 1000, shape = shape, scale = scale)
plot(density(invg))
for(i in 1:20){
  abline(v=1/(rgamma(1, shape=shape, scale=scale)))
}


#InvG scale = scale 
#G scale = rate, 1/()
#result - better
shape=5
scale=4
rate=1/scale
invg <- rinvgamma(n = 1000, shape = shape, scale = scale)
plot(density(invg))
for(i in 1:20){
  abline(v=1/(rgamma(1, shape=shape, scale=rate)))
  abline(v=1/(randGamma(shape, 1/rate)), col="red")
}




#TEST Inverse Gamma
pdf("invGammaTest.pdf")
shape <- 8
rate <- 2
x <- rinvgamma(n=1000, shape=shape, scale=rate)
plot(density(x))
for(i in 1:20)
{
  abline(v = rinvgamma(1, shape=shape, scale=rate))
  abline(v=sqrt(1/(randGamma(shape, rate))), col="red")
  abline(v=sqrt(1/ rgamma(1, shape=shape, scale=rate)), col="green")
}
dev.off()

#TEST randNorm
pdf("randNormTest.pdf")
mean <- 3
sd <- 1
x <- rnorm(1000, mean = mean, sd = sd)
plot(density(x))
for(i in 1:20)
{
  abline(v=randNorm(mean, sd))
}
dev.off()


#TEST randLogNorm
pdf("randLogNormTest.pdf")
mean <- 3
sd <- 1
x <- rlnorm(1000, mean = mean, sd = sd)
plot(density(x))
for(i in 1:20)
{
  abline(v=randLogNorm(mean, sd))
}
dev.off()


#TEST randExp
pdf("randExpTest.pdf")
rate <- 2
x <- rexp(1000, rate = rate)
plot(density(x))
for(i in 1:20)
{
  abline(v=randExp(rate))
}
dev.off()



#TEST randUnif
pdf("randUnifTest.pdf")
min <- 3
max <- 17
x <- runif(1000, min = min, max = max)
plot(density(x))
for(i in 1:20)
{
  abline(v=randUnif(mix, max))
}
dev.off()


