# IntroLikelihoodToBayes

![](https://github.com/ManueleBazzichetto/IntroLikelihoodToBayes/blob/main/Bayes_pic.jpeg)

This repo includes slides and an R script (in the Examples folder) used on the 1st day of the course "Introduction to Likelihood and Bayes" (by Manuele Bazzichetto and Petr Keil) held at the Department of Spatial Sciences of the Czech University of Life Sciences.

The material covers topics that range from data to parametric models for statistical populations (DataModels presentation), and from basic principles of probability and maximum likelihood estimation to Bayes' rule and Bayesian estimation of model parameters (ProbabilityLikelihoodBayes).

Slides are provided as both html and .qmd files, so that people will be able to modify, correct, integrate any material presented during the course. 

For most of the R examples (e.g., LikelihoodEstimation.R), I used the (amazing) palmerpenguins dataset[^1], which is made available through a namesake R package.

The "Examples" folder includes 4 practical examples:

- brute force maximum likelihood estimation of the population mean of a Gaussian model;
- same task as above using the optim() R function (and the "BFGS" method);
- maximum likelihood estimation of the population rate of a Poisson model;
- approximation grid bayesian estimation of the population rate of a Poisson model.

For any doubt or inquiry, please contact me at: manuele.bazzichetto@gmail.com

[^1]: Horst AM, Hill AP, Gorman KB (2020). palmerpenguins: Palmer Archipelago (Antarctica) penguin data. R package version 0.1.0. https://allisonhorst.github.io/palmerpenguins/. doi: 10.5281/zenodo.3960218.


