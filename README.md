# piecewiseSEM: Piecewise Structural Equation Modeling in R
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/piecewiseSEM)](https://cran.r-project.org/package=piecewiseSEM)

## Version 2.3.01
## Last updated: 01 June 2023

## To install

Run the following code to install the latest version from CRAN:
```
install.packages("piecewiseSEM")
```
Run the following code to install the development version:
```
devtools::install_github("jslefche/piecewiseSEM@devel")
```
Note: the development version may be unstable and lead to unanticipated bugs.
Contact the package developer with any bugs or issues.

## Getting Help
See our website at [piecewiseSEM](http://jslefche.github.io/piecewiseSEM/)

There is an online resource available for SEM, including `piecewiseSEM` and `lavaan`, available [https://jslefche.github.io/sem_book/](https://jslefche.github.io/sem_book/)

Version 2 is a major update to the `piecewiseSEM` package that uses a completely revised syntax that better reproduces the base R syntax and output. It is highly recommended that consult the resource above even if you have used the package before as it documents the many changes.

Currently supported model classes: `lm, glm, gls, Sarlm, lme, glmmPQL, lmerMod, merModLmerTest, glmerMod. glmmTMB, gam`

### Example
```
# Load library
library(piecewiseSEM)

# Create fake data
set.seed(1)

data <- data.frame(
  x = runif(100),
  y1 = runif(100),
  y2 = rpois(100, 1),
  y3 = runif(100)
)

# Create SEM using `psem`
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
