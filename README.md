# piecewiseSEM: Piecewise Structural Equation Modeling in R

## Version 2.0

This is a major update to the `piecewiseSEM` package that includes new functionality, and a completely revised syntax that better reproduces the base R syntax and output.

Currently supported model classes: `lm, glm, gls, lme, glmmPQL, lmerMod, merModLmerTest, glmerMod`

###WARNING: THIS IS A HIGHLY EXPERIMENTAL BRANCH. FOR DEVELOPMENT ONLY: DO NOT USE!!

### Example
```
# Install from github
library(devtools)
install_github("jslefche/piecewiseSEM@2.0")

# Load library
library(piecewiseSEM)

# Create fake data
set.seed(1) 

data = data.frame(
  x = runif(100),
  y1 = runif(100),
  y2 = rpois(100, 1),
  y3 = runif(100)
)

# Store in SEM list 
modelList = psem(
  lm(y1 ~ x, data),
  glm(y2 ~ x, "poisson", data),
  lm(y3 ~ y1 + y2, data)
)

# Run summary
summary(modelList, data)

# Address conflict using conserve = T
summary(modelList, data, conserve = T)

# Address conflict using direction = c()
summary(modelList, data, direction = c("y2 <- y1"))

# Address conflict using correlated errors
modelList2 = update(modelList, y2 %~~% y1)

summary(modelList2, data)
```
