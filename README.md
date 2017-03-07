# piecewiseSEM: Piecewise Structural Equation Modeling in R

## Version 2.0

###WARNING: THIS IS A HIGHLY EXPERIMENTAL BRANCH. FOR DEVELOPMENT ONLY: DO NOT USE!!

### Example
```
# Install from github
library(devtools)
install_github("jslefche/piecewiseSEM@2.0")

# Load library
library(piecewiseSEM)

# Create fake data
data = data.frame(
  x = runif(100),
  y1 = runif(100),
  y2 = rpois(100, 1),
  y3 = runif(100)
)

# Store in SEM list 
modelList = list.sem(
  lm(y1 ~ x, data),
  glm(y2 ~ x, "poisson", data),
  lm(y3 ~ y1 + y2, data)
)

# Run summary
summary(modelList)

# Address conflict using conserve = T
summary(modelList, conserve = T)

# Address conflict using direction = c()
summary(modelList, direction = c("y2 <- y1"))
```
