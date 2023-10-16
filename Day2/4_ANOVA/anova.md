---
title: "ANOVA -- **part 1**: The data and fixed-effects model definition"
author: "Petr Keil"
date: "October 2023"
output:
  html_document:
    highlight: pygments
    keep_md: yes
    number_sections: yes
    theme: cerulean
    toc: yes
  pdf_document: default
---

***

# Objective

The aim of this lesson is to leave the participants to come up with their
code for simple one-way ANOVA (part 1), and to experiment with random effects ANOVA (part 2).

***

# The Data

We will use modified data from the example from **Marc Kery's Introduction to WinBUGS for Ecologists**, page 119 (Chapter 9 - ANOVA). The data describe snout-vent lengths in 5 populations of Smooth snake (*Coronella austriaca*).

![](figure/snake.png)

***

Loading the data from the web:


```r
  snakes <- read.csv("https://raw.githubusercontent.com/petrkeil/ML_and_Bayes_2021_CZU/main/09_ANOVA/snakes.csv")

  summary(snakes)
```

```
##    population   snout.vent   
##  Min.   :1    Min.   :33.19  
##  1st Qu.:2    1st Qu.:45.11  
##  Median :3    Median :48.56  
##  Mean   :3    Mean   :50.45  
##  3rd Qu.:4    3rd Qu.:58.09  
##  Max.   :5    Max.   :65.49
```

Plotting the data:

```r
  par(mfrow=c(1,2))
  plot(snout.vent ~ population, data=snakes,
       ylab="Snout-vent length [cm]")
  boxplot(snout.vent ~ population, data=snakes,
          ylab="Snout-vent length [cm]",
          xlab="population",
          col="grey")
```

![](anova_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

***

# The model: Fixed-effects ANOVA

For a given snake $i$ in population $j$ **the model** can be written as:

$y_{ij} \sim Normal(\alpha_j, \sigma)$

***

# Tasks for you:

* Write this model in the BUGS language and dump it into a file
using `cat`.

* Pepare the data for this model to the `list` format.

* Fit the model and estimate posterior distributions of $\alpha_j$.

* Is there a significant difference of mean snout-vent length between 
populations 1 and 2?





