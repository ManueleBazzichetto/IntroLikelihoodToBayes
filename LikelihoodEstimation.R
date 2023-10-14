##Maximum likelihood estimation

library(ggplot2)
library(ggpubr)
library(palmerpenguins)
library(manipulate)

#To install missing packages use install.packages()

#import penguins dataset
data("penguins")

#coerce to data.frame (you don't have to, it's me that I am old-fashion)
penguins <- as.data.frame(penguins)

#have a look at what's inside penguins
str(penguins)
head(penguins)

#subset penguins to only keep data on Gentoo
gentoo_df <- penguins[penguins$species == "Gentoo", ]

#we have 1 NA in the column for the body mass, which is the measure we are interested in
sum(is.na(gentoo_df$body_mass_g)) #1

#let's get rid of it
gentoo_df <- gentoo_df[!is.na(gentoo_df$body_mass_g), ]

anyNA(gentoo_df$body_mass_g) #FALSE

gentoo_body_mass <- gentoo_df$body_mass_g

##------------------1st example!
##------------------estimate population mean of Gentoo body mass (BM) through maximum likelihood estimation

#first of all, we will assume (for the sake of simplicity) that we know the population variance of Gentoo body mass
#this will reduce the estimation problem to a 1-dimension

#we will use a 'brute force' approach to estimate the population mean of BM
#this basically means that we will generate a sequence of "potential" candidates for the BM population mean
#we will then derive the joint probability of our sample of BM given each of possible candidate for the BM population mean
#once we have the likelihood function, we pick up the candidate associated with the max value of the joint probability of the data
#this will be the maximum likelihood estimator

#create a sequence of candidate values for the BM population mean

#let's find out what's the sample mean
BM_sample_mean <- mean(gentoo_body_mass) #5076.016 gr

#we can set a sequence of possible means that vary from 4500 (gr) to 5500 (gr) by increments of 10 gr
BM_mean_seq <- seq(from = 4500, to = 5500, by = 10)

#let's now assume we know the population variance, which we fix at the sample variance
#we actually compute the sqrt of the variance, i.e. the standard deviation, as we will need that later
BM_sample_sd <- sd(gentoo_body_mass)

#how do we compute the likelihood?
#we will use the dnorm function, which, when fed with some values for the mean and sd parameters of the Gaussian, will return the 
#corresponding density for each entry of the sample
#summing up all densities will give us the (log)likelihood associated with the given mean and sd
#setting the log argument equal to TRUE will return the log-transformed likelihood, as the density of each data point will be log-transformed 


#let's now write a function which will compute the joint probability of our sample of BM given the potential values of the mean
BruteForce_likelihood <- sapply(BM_mean_seq, function(i) {
  LL_value <- sum(dnorm(gentoo_body_mass, mean = i, sd = BM_sample_sd, log = TRUE))
  return(LL_value)
})

#let's see what we obtained
BruteForce_df <- data.frame(LLikelihood = BruteForce_likelihood, Mean_values = BM_mean_seq)

ggplot(BruteForce_df, aes(x = Mean_values, y = LLikelihood)) +
  geom_point() +
  theme_pubr() 

#so, who's the winner? Let's find the maximum likelihood estimator, i.e. the potential candidate for the mean 
#that maximized the log-likelihood
BruteForce_df[which.max(BruteForce_df$LLikelihood), "Mean_values"] #5080

max_ll_value <- BruteForce_df$LLikelihood[which.max(BruteForce_df$LLikelihood)]

#let's zoom in and see the winner
ggplot(BruteForce_df, aes(x = Mean_values, y = LLikelihood)) +
  geom_point() +
  geom_point(data = BruteForce_df[which.max(BruteForce_df$LLikelihood), ],
             aes(x = Mean_values, y = LLikelihood), col = "red", size = 4) +
  ylim(max_ll_value - 2, max_ll_value) +
  xlab("Candidates for the mean") +
  theme_pubr() +
  theme(text = element_text(size = 18))


#so the MLE is 5080!
#note that this value is very close to the sample mean!
#funny enough, this is because the maximum likelihood estimator of the population mean of a Gaussian is the sample mean!
#so, why didn't we get the sample mean as the maximum likelihood estimator?
#this is due to the 'brute force' approach, which forces us to set a sequence of possible mean values and also a resolution
#of the sequence (the incremental steps). So, simply put, the sample mean was not among the potential candidates for the mean


##------------------same exercise, but using an interactive interface

manipulate(
  {
    LL_val <- sum(dnorm(gentoo_body_mass, mean = MeanCandidate, sd = BM_sample_sd, log = TRUE))
    par(mfrow = c(1, 2))
    Gauss_densities <- dnorm(seq(min(gentoo_body_mass), max(gentoo_body_mass), by = 1),
                             mean = MeanCandidate, sd = BM_sample_sd)
    hist(gentoo_body_mass, xlab = "Gentoo's body mass - sample", 
         breaks = 15, main = NULL, freq = F)
    lines(seq(min(gentoo_body_mass), max(gentoo_body_mass), by = 1),
          Gauss_densities, col = "black", lwd = 2)
    abline(v = MeanCandidate, col = "green", lwd = 2)
    plot(BM_mean_seq, BruteForce_df$LLikelihood,
         xlab = "Candidates for the mean", ylab = "LL")
    points(MeanCandidate, LL_val, col = "red", pch = 16, cex = 2)
    abline(v = MeanCandidate, col = "red", lwd = 2)
    text(x = 5400, y = -940, labels = paste("LL = ", round(LL_val, 2)))
    },
  
  MeanCandidate = slider(min = min(BM_mean_seq), max = max(BM_mean_seq), initial = BM_sample_mean,
                   label = "Candidate for the mean")
)

dev.off()

#there are algorithms that find the MLE for us!
#this comes especially handy when we need to maximize more complex likelihood functions
#(e.g., multiple parameters or hierarchical models)


##------------------2nd example!
##---------------estimate regression parameters for the relationship between Gentoo body mass (BM) and flipper length (FL)
##---------------through maximum likelihood estimation

#quickly check if there are NAs in the column for FL
anyNA(gentoo_df$flipper_length_mm) #FALSE

#the scientific question itself is not interesting, but it can give us an idea of how we estimate regression parameters using the likelihood

#first of all, let's have a look at the relationship between BM and flipper length (FL) in the sample
#as usual, this simply gives us an idea of the possible relationship between these two measures in the population
#as usual, we should have a (possibly parametric) model for the population in mind

#note that we put body_mass_g (BM) on the y-axis, and we will use it as the response variable in the regression, 
#but this does not mean that flipper_length_mm (FL), which we will use as the predictor, actually 'causes' body_mass_g
#we could do the opposite and see the association between these two measures does not change (see correlation)
ggplot(gentoo_df, aes(y = body_mass_g, x = flipper_length_mm)) +
  geom_point(alpha = .6) +
  ylab("Body mass (BM)") + xlab("Flipper length (FL)") +
  theme_pubr()

#no surprise that BM increases as FL increases! We knew this was not a stimulating research question..

#let's anyway dive into the regression problem..

#What's regression about - long explanation
#####---------- 

#as we saw during the 'theoretical' part of this class, regressing BM on FL is not only about finding a line that goes 
#across the data-points and gives us an idea of the association between two measures

#what we actually do is to fit up to N (N = sample size) parametric models, assuming that (at each combination of response
#and predictor(s) values) the scatter of the measures around their (population) mean behaves in a particular way
#specifically, we assume that: 
#1) this scatter is centered at 0, basically meaning the scatter is white noise around the mean;
#2) this scatter is homogeneous 'along' the regression line;
#3) (if not assumed differently) the scatter at a given combination of response/predictor value does not depend on 
#the scatter at another combination

#if we further assume this scatter is Gaussian around the population mean, we have a parametric model to use to maximize the likelhood!

#we can visualize this situation as up to N Gaussian distributions popping up along the regression line, which gives us 
#the value of the population mean for each of these Gaussian distributions
#(if you can't visualize this, go back to the 'DataModels' slides)

#if we have Gaussian(s) distributions, this means we have population mean(s) and standard deviation(s) to estimate!
#however, this time the mean will depend on the value of FL, as we are assuming we can predict BM given the value of FL

#basically: mean BMi = alpha + beta * FLi
#where alpha is the intercept (mean BMi when FLi is 0, pretty unlikely event) of the regression line, 
#and beta is the slope of the regression line (which literally gives us the mean change in BM for a 1-unit increment in FLi)

#what about the standard deviation(s)? Well, we assume the scatter around the mean is homogeneous,
#so let's assume the sd(s) are all the same

#regression parameters are usually estimated using ordinary least squares (OLS), but, given we have a parametric model for the
#behavior of our measures along the mean, we can use maximum likelihood estimation
#(note that, under the Gaussian assumptions, OLS and maximum likelihood estimation lead to the same answer,
#except for the estimation of the sd)

####-----------


#to estimate the regression parameters for the relationship between BM and FL, as well as the standard deviation(s) 
#measuring the scatter of BM around the population mean, we will not use brute force, while rather an algorithm that will 
#climb the likelihood surface for us. It will find the parameter values that maximize the joint probability
#of our sample given the model

#this can be achieved using the optim() R function.

#first of all, we need to write a short R function which takes BM, FL and a vector of parameters as input
#and returns the log-likelihood evaluated at each combination of parameter values

RegressionLikelihood_estim <- function(x, BM, FL_cent) {
  #x is an hypothetical vector with length equal to the number of params we need to estimate
  #we subset the elements of this hypothetical vector and assign them to specific objects:
  #intercept of the regression line
  alpha <- x[1]
  #slope of the regression line
  beta <- x[2]
  #scatter of the measures around the mean
  std_dev <- x[3]
  #log-likelihood function - note that we take the negative of the log-likelihood as optim() minimizes functions
  #over the parameter space
  ll <- -sum(dnorm(BM, mean = alpha + beta * FL_cent, sd = std_dev, log = T))
  return(ll)
}

#we need to set an initial value for each of the parameters (passed as the first argument in optim) 
#this will help the algorithm converge fast

##--initial parameter values

#beta: let's say that on average BM increases of 50 gr for each 1-mm increment in FL
#alpha: we center FL, so that alpha is the mean value of BM when FL is also at its mean value

#center FL
FL_cent <- with(gentoo_df, flipper_length_mm - mean(flipper_length_mm))
#sanity check
mean(FL_cent)

#mean value of FL
mean(gentoo_df$flipper_length_mm) #217.187
#mean of some BM values around the mean value of FL
mean(gentoo_df[gentoo_df$flipper_length_mm < 218 & gentoo_df$flipper_length_mm > 216, "body_mass_g"]) #5016.667

#sd: let's say that the average scatter of BM around the regression line is of about 200 gr

##--

#run the optim function! We use the "BFGS" method, which uses derivatives to minimize the negative log-likelihood
RegressionLikelihood_output <- optim(par = c(5016, 50, 200), fn = RegressionLikelihood_estim, method = "BFGS",
                                     BM = gentoo_body_mass, FL_cent = FL_cent, hessian = T)


#let's scan the output

#estimated parameters - aka maximum likelihood estimates
BM_FL_params <- RegressionLikelihood_output$par #our initial guess about the value of the parameters was not that bad

BM_FL_params <- setNames(BM_FL_params, c("int", "slope", "sd"))

#value of the log-likelihood evaluated at the parameters value above
RegressionLikelihood_output$value

#diagnostic - if the value of convergence is 0, then the algorithm successfully converged to a solution
RegressionLikelihood_output$convergence

#we also got the Hessian matrix for free, which is a matrix of second order partial derivatives of the log-likelihood 
#function evaluated around the maximum likelihood estimates

#taking the inverse of this function will give us a measure of the precision of the maximum likelihood estimators
#basically, an (asymptotically valid) covariance matrix of the maximum likelihood estimators
solve(RegressionLikelihood_output$hessian)


#let's compare what we got using optim with what we obtain using a more 'conventional' tool to estimate
#regression parameters under the Gaussian model

coef(lm(gentoo_body_mass ~ FL_cent)) #intercept: 5076.0163; slope: 54.6225
sigma(lm(gentoo_body_mass ~ FL_cent)) #360.1676

coef(glm(gentoo_body_mass ~ FL_cent, family = "gaussian")) #intercept: 5076.0163; slope: 54.6225
sqrt(15696203/length(gentoo_body_mass)) #357.2274

#if we use less data, we see that we get pretty different results for the standard deviation(s)
#this is because the maximum likelihood estimator of the variance is slightly downwardly biased
indexes_to_sample <- sample(length(gentoo_body_mass), 5, replace = F)

sigma(lm(gentoo_body_mass[indexes_to_sample] ~ FL_cent[indexes_to_sample])) #223.8

summary(glm(gentoo_body_mass[indexes_to_sample] ~ FL_cent[indexes_to_sample], family = "gaussian")) 
sqrt(150279/5) #173.3661


#now that we have our parameters, we can try to see the estimated Gaussian(s) along the regression line

#we'll put all quantities we are interested in into a data.frame, so that we can easily plot them

BM_FL_result <- data.frame(BM = gentoo_body_mass, FL = gentoo_df$flipper_length_mm,
                           EstMean = BM_FL_params[["int"]] + BM_FL_params[["slope"]] * FL_cent,
                           EstScatter_low = BM_FL_params[["int"]] + BM_FL_params[["slope"]] * FL_cent - 2*BM_FL_params[["sd"]],
                           EstScatter_up = BM_FL_params[["int"]] + BM_FL_params[["slope"]] * FL_cent + 2*BM_FL_params[["sd"]])


ggplot(BM_FL_result, aes(x = FL, y = BM)) +
  geom_ribbon(data = BM_FL_result, aes(ymin = EstScatter_low, ymax = EstScatter_up), alpha = .3) +
  geom_point(alpha = .6) +
  geom_line(data = BM_FL_result, aes(x = FL, y = EstMean), col = "green", lwd = 1.2) +
  ylab("Body mass") + xlab("Flipper length") +
  theme_pubclean() +
  theme(text = element_text(size = 18))


##------------------3rd example!
##---------------estimate regression parameters for the relationship between Gentoo puppies and adult body mass (BM)
##---------------through maximum likelihood estimation


#we don't have information on the number of puppies per Gentoo, but we can simulate this data

#we will do as follows:
#first, we define a model for the number of puppies a Gentoo adult generates during her/his life
#we assume the number of puppies is Poisson distributed with rate parameter lambda
#we further assume that the number of of puppies increases as a function of the adult's body mass
#basically, the 'bigger' the adult, the higher the number of puppies
#also, as predicted by the Poisson model, the variability in the number of puppies increases with the mean
#basically, 'bigger' penguins can make lots or few puppies

#here, for the sake of simplicity (and the lack of data), we equally treat males and females
#obviously, this assumption does not make any sense

#let's convert BM from grams to kg
BM_cent <- gentoo_body_mass/1000

#then, let's center BM expressed in kg
BM_cent <- BM_cent - mean(BM_cent)

#the model for the rate parameter is lambda = exp(alpha + beta * BM_cent)
#we need the exp() to be sure we don't predict negative number of puppies for adults with low body mass

##--parameters

alpha_pois <- 1.7

beta_pois <- 0.45

#we set alpha = 1.7, which, given that we centered BM, means that on the counts scale the average number of puppies generated
#by an adult of average body mass is exp(1.7) = 5.5

#we set beta = 0.45 on the link scale, which means that on the counts scale, a 1-unit kg increase in the adult 
#translates in a exp(0.45) = 56% increase of the average number of puppies

#let's see the predicted range of lambda parameters:
lambda_pois <- exp(alpha_pois + beta_pois * BM_cent)

range(lambda_pois)

##--

#we can now simulate the number of puppies generated by each adult given her/his body mass

NumberPuppies_df <- data.frame(NumPuppies = rpois(length(BM_cent), lambda = lambda_pois),
                               BM = BM_cent,
                               TrueLambdas = lambda_pois)


#let's see what we got
ggplot(NumberPuppies_df, aes(x = BM, y = NumPuppies)) +
  geom_point(alpha = 0.6) +
  geom_line(aes(x = BM, y = TrueLambdas), col = "green", lwd = 1.2) +
  ylab("Number of puppies") + xlab("Body mass centered") +
  theme_pubclean() 


#we can now use optim() to estimate the regression parameters using maximum likelihood estimation

PoissonLikelihood_estim <- function(x, NP, BM_cent) {
  alpha <- x[1]
  beta <- x[2]
  #log-likelihood function - note that we take the negative of the log-likelihood as optim() minimizes functions
  #over the parameter space
  ll <- -sum(dpois(NP, lambda = exp(alpha + beta * BM_cent), log = T))
  return(ll)
}

#we need to provide initial parameter values. This time we know them..
#let's anyway make it a bit hard for the algorithm. We'll provide values that slightly differ from the true ones
PoissonLikelihood_output <- optim(par = c(1.5, 0.40), fn = PoissonLikelihood_estim, method = "BFGS",
                                     NP = NumberPuppies_df$NumPuppies,
                                  BM_cent = NumberPuppies_df$BM, hessian = T)


#we can now compare the number of puppies estimated by the model vs. the simulated ones
NumberPuppies_df$EstLambdas <- with(NumberPuppies_df, exp(PoissonLikelihood_output$par[1] + PoissonLikelihood_output$par[2] * BM))

ggplot(NumberPuppies_df, aes(x = BM, y = NumPuppies)) +
  geom_point(alpha = 0.6) +
  geom_line(aes(x = BM, y = TrueLambdas), col = "green", lwd = 1.2) +
  geom_line(aes(x = BM, y = EstLambdas), col = "red", lwd = 1.2) +
  ylab("Number of puppies") + xlab("Body mass centered") +
  theme_pubclean() 

#we also got the Hessian matrix for free, which is a matrix of second order partial derivatives of the log-likelihood 
#function evaluated around the maximum likelihood estimates

#taking the inverse of this function will give us a measure of the precision of the maximum likelihood estimators
#basically, an (asymptotically valid) covariance matrix of the maximum likelihood estimators
solve(PoissonLikelihood_output$hessian)

vcov(glm(NumPuppies ~ BM, data = NumberPuppies_df, family = 'poisson'))

