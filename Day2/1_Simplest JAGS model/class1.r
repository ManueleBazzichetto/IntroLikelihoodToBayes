snakes <- read.csv("https://raw.githubusercontent.com/petrkeil/ML_and_Bayes_2021_CZU/main/09_ANOVA/snakes.csv")

summary(snakes)

par(mfrow=c(1,2))
plot(snout.vent ~ population, data=snakes,
     ylab="Snout-vent length [cm]")
boxplot(snout.vent ~ population, data=snakes,
        ylab="Snout-vent length [cm]",
        xlab="population",
        col="grey")

snake.data <- list(y=snakes$snout.vent,
                   x=snakes$population,
                   N=nrow(snakes), 
                   N.pop=5)


cat("

model
{
    # priors
    grand.mean ~ dnorm(0, 0.001)
    grand.sigma ~ dunif(0,100)
    grand.tau <- 1/pow(grand.sigma)
  
    sigma ~ dunif(0,100)
    tau <- 1/(sigma*sigma)
  
    for(j in 1:N.pop)
    {
      alpha[j] ~ dnorm(grand.mean, grand.tau)
    }
    
    # likelihood
    for(i in 1:N)
    {
      y[i] ~ dnorm(alpha[x[i]], tau)
    }
  
    # derived quantity
    delta12 <- alpha[1] - alpha[2]
  
}    
    
", file = "anova.txt")




model.fit.fix <- jags(data=snake.data, 
                      model.file="anova.txt",
                      parameters.to.save=c("alpha", "delta12"),
                      n.chains=3,
                      n.iter=2000,
                      n.burnin=1000,
                      DIC=FALSE)

library(mcmcplots)
caterplot(model.fit.fix)
