# IntroLikelihoodToBayes

![](https://github.com/ManueleBazzichetto/IntroLikelihoodToBayes/blob/main/Bayes_pic.jpeg)

[Thomas Bayes](https://en.wikipedia.org/wiki/Thomas_Bayes)

This repo includes slides and an R script (in the Examples and Day2 folder) used for the 2-day course on "Introduction to Likelihood and Bayes" (by Manuele Bazzichetto and Petr Keil) held at the Department of Spatial Sciences of the Czech University of Life Sciences.

Day1: The material covers topics that range from data to parametric models for statistical populations (DataModels slides), and from basic principles of probability and maximum likelihood estimation to Bayes' rule and Bayesian estimation of model parameters (ProbabilityLikelihoodBayes slides). Slides are provided as both html and .qmd files, so that people will be able to modify, correct, integrate any material presented during the course. 

For most of the Day1 R examples (e.g., LikelihoodEstimation.R), we used the (amazing) palmerpenguins dataset[^1], which is made available through a namesake R package.

The "Examples" folder includes 4 practical examples:

- brute force maximum likelihood estimation of the population mean of a Gaussian model;
- using the optim() R function (and the "BFGS" method) to estimate parameters of a Gaussian regression;
- using the optim() R function (and the "BFGS" method) to estimate parameters of a Poisson regression;
- approximation grid bayesian estimation of the population rate of a Poisson model.

Day2: Practical examples of Bayesian estimation and inference using JAGS.

For any doubt or inquiry, please contact me at: manuele.bazzichetto@gmail.com

[^1]: Horst AM, Hill AP, Gorman KB (2020). palmerpenguins: Palmer Archipelago (Antarctica) penguin data. R package version 0.1.0. https://allisonhorst.github.io/palmerpenguins/. doi: 10.5281/zenodo.3960218.


