###### Build script
library(devtools)
setwd("../")
# library(Rd2roxygen)
# Rd2roxygen::Rd2roxygen("piecewiseSEM")


piecewiseSEM <- as.package("./piecewiseSEM") 
devtools::use_build_ignore(c("build.R", ".git", ".gitignore", "R/anova.R"),
                           pkg = "./piecewiseSEM")

### Build documentation and vignettes
document(piecewiseSEM)
#clean_vignettes(piecewiseSEM)
#build_vignettes(piecewiseSEM)


### Load package
load_all(piecewiseSEM, reset=T)

### Check and build
#check(piecewiseSEM, cran=TRUE)
#build(piecewiseSEM, path="./piecewiseSEM/builds")
#devtools::check_built("./piecewiseSEM_2.0.tar.gz")
install(piecewiseSEM) 
#run_examples(multifunc) 

### Build pkgdown site
library(pkgdown)
setwd("piecewiseSEM")
build_site()



### Examples to verify
library(piecewiseSEM)
data(keeley)

mod <- psem(
  lm(rich ~ cover, data=keeley),
  lm(cover ~ firesev, data=keeley),
  lm(firesev ~ age, data=keeley),
  data = keeley
  
)

res <- residuals(mod)
anova(mod)
fisherC(mod)

mod2 <- psem(
  lm(rich ~ cover, data=keeley),
  lm(cover ~ firesev + age, data=keeley),
  lm(firesev ~ age, data=keeley),
  rich %~~% firesev,
  data = keeley
  
)

anova(mod, mod2)


#test mixed models
library(lme4)
library(nlme)

# Create list of structural equations
shipley_psem <- psem(

  lme(DD ~ lat, random = ~ 1 | site / tree, na.action = na.omit,
  data = shipley),

  lme(Date ~ DD, random = ~ 1 | site / tree, na.action = na.omit,
  data = shipley),

  lme(Growth ~ Date, random = ~ 1 | site / tree, na.action = na.omit,
  data = shipley),

  glmer(Live ~ Growth + (1 | site) + (1 | tree),
  family = binomial(link = "logit"), data = shipley)

  )

summary(shipley_psem)



# Create list of structural equations
shipley_psem_lme4 <- psem(
  
  lmer(DD ~ lat + (1 | site / tree), 
      data = shipley),
  
  lmer(Date ~ DD + (1 | site / tree), 
      data = shipley),
  
  lmer(Growth ~ Date + (1 | site / tree),
      data = shipley),
  
  glmer(Live ~ Growth + (1 | site) + (1 | tree),
        family = binomial(link = "logit"), data = shipley),
  
  data = shipley
  
)

summary(shipley_psem_lme4)

#lmerTest
library(lmerTest)
shipley_psem_lmerTest <- psem(
  
  lmer(DD ~ lat + (1 | site / tree), 
       data = shipley),
  
  lmer(Date ~ DD + (1 | site / tree), 
       data = shipley),
  
  lmer(Growth ~ Date + (1 | site / tree),
       data = shipley),
  
  glmer(Live ~ Growth + (1 | site) + (1 | tree),
        family = binomial(link = "logit"), data = shipley),
  
  data = shipley
  
)

summary(shipley_psem_lmerTest)
