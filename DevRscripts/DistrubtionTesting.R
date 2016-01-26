library(ribModel)



#TEST randGamma
pdf("randGammaTest.pdf")
shape <- 1
rate <- 2
x <- rgamma(1000, shape = shape, rate = rate)
plot(density(x))
for(i in 1:20)
{
  abline(v=randGamma(shape, rate))
}
dev.off()


#TEST Inverse Gamma
pdf("invGammaTest.pdf")
shape <- 100
rate <- 20
x <- rinvgamma(n=1000, shape=shape, scale=1/rate)
plot(density(x))
for(i in 1:20)
{
  abline(v=sqrt(1/(randGamma(shape, rate))))
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


