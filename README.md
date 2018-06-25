# piecewiseSEM: Piecewise Structural Equation Modeling in R

## Version 2.0.1
## Last updated: 25 June 2018

This version is a major update to the `piecewiseSEM` package that usesa completely revised syntax that better reproduces the base R syntax and output. It is highly recommended that consult `vignette("piecewiseSEM")` even if you have used the package before as it documents the many changes.

It also incorporates new functionality in the form of coefficient standardization and updated methods for R^2 for mixed models. 

Currently supported model classes: `lm, glm, gls, pgls, sarlm, lme, glmmPQL, lmerMod, merModLmerTest, glmerMod`

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
