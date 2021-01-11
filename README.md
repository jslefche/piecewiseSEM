# piecewiseSEM: Piecewise Structural Equation Modeling in R
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/piecewiseSEM)](https://cran.r-project.org/package=piecewiseSEM)

## Version 2.2.0
## Last updated: 11 January 2021

## To install

Run the following code to install the development version:
```
# install.packages("devtools") # if devtools not installed

devtools::install_github("jslefche/piecewiseSEM@devel")
```

## Getting Help
See our website at (piecewiseSEM)[http://jslefche.github.io/piecewiseSEM/]

There is an online resource available for SEM, including `piecewiseSEM` and `lavaan`, available (here)[https://jslefche.github.io/sem_book/]

This version is a major update to the `piecewiseSEM` package that usesa completely revised syntax that better reproduces the base R syntax and output. It is highly recommended that consult `vignette("piecewiseSEM")` even if you have used the package before as it documents the many changes.

It also incorporates new functionality in the form of coefficient standardization and updated methods for R^2 for mixed models.

Currently supported model classes: `lm, glm, gls, lme, glmmPQL, lmerMod, merModLmerTest, glmerMod`

## To add a new model class
Frequently, we get requests to add new model classes. We'd like to accomodate wherever we can! Currently, to add a new model class, you will need to update the following functions in the following files:  
- In `coefs.R`
      - `getCoefficients()` - will need a method to generate a standardized coefficient table from your class
      - `GetSDy()` - will need a method to get the proper SD of y, particularly for non-Gaussian error families
- In `helpers.R`
      - `GetSingleData()` - need a method to get a data frame from your class type
      - `GetOLRE()` - if your class has random effects, will need a method to get observation level random effects
      - `nObs()` - if your class has a non-standard way to get the number of observations, will need it here.
      - You might need to write a `all.vars` method for you class. See `all.vars.merMod()` as an example.
- In `psem.R`
      - `evaluateClasses()` - so we know that this class can work in `piecewiseSEM`
- In `rsquared.R`
      - `rsquared()` - what function should be used to get the rsquared for your model class?
      - Will need an `rsquared.youclass()` method to get the R^2

### Example
```
# Install development branch from github
library(devtools)
install_github("jslefche/piecewiseSEM@devel", build_vignette = TRUE)

# Load library
library(piecewiseSEM)

# Read vignette
vignette("piecewiseSEM")

# Create fake data
set.seed(1)

data <- data.frame(
  x = runif(100),
  y1 = runif(100),
  y2 = rpois(100, 1),
  y3 = runif(100)
)

# Store in SEM list
modelList <- psem(
  lm(y1 ~ x, data),
  glm(y2 ~ x, "poisson", data),
  lm(y3 ~ y1 + y2, data),
  data
)

# Run summary
summary(modelList)

# Address conflict using conserve = T
summary(modelList, conserve = T)

# Address conflict using direction = c()
summary(modelList, direction = c("y2 <- y1"))

# Address conflict using correlated errors
modelList2 <- update(modelList, y2 %~~% y1)

summary(modelList2)
```
