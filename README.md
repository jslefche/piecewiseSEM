# piecewiseSEM: Piecewise Structural Equation Modeling in R

## Version 2.0
## Last updated: 24 October 2017

This is a major update to the `piecewiseSEM` package that includes new functionality, and a completely revised syntax that better reproduces the base R syntax and output. It is advised that consult `vignette("piecewiseSEM")` even if you have used the package before as it documents the many changes.

Currently supported model classes: `lm, glm, gls, pgls, sarlm, lme, glmmPQL, lmerMod, merModLmerTest, glmerMod`

### Example
```
# Install from github
library(devtools)
install_github("jslefche/piecewiseSEM@2.0")

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
