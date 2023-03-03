#####-----Build script

# generate documentation
roxygen2::roxygenise()

# build documentation and vignettes
document(piecewiseSEM)
#clean_vignettes(piecewiseSEM)
#build_vignettes(piecewiseSEM)

# create website for package
pkgdown::build_site()

# build package
piecewiseSEM <- as.package("./piecewiseSEM") 
devtools::use_build_ignore(c("build.R", ".git", ".gitignore", "docs"),
                           pkg = "./piecewiseSEM")
# Load package
load_all(piecewiseSEM, reset = T)

### Check and build
check(piecewiseSEM, cran=TRUE)
build(piecewiseSEM, path="./piecewiseSEM/builds")
devtools::check_built("./piecewiseSEM_2.0.tar.gz")
install(piecewiseSEM) 
#run_examples(piecewiseSEM) 

### Build pkgdown site
library(pkgdown)
setwd("piecewiseSEM")
build_site()


#####-----Examples to verify functionality
library(piecewiseSEM)
data(keeley)

# fit model
mod <- psem(
  lm(rich ~ cover, data=keeley),
  lm(cover ~ firesev, data=keeley),
  lm(firesev ~ age, data=keeley),
  data = keeley)

mod

# d-sep tests
basisSet(mod)
dSep(mod)
fisherC(mod)
AIC(mod, AIC.type = "dsep")

# Chisq
LLchisq(mod)
AIC(mode, AIC.type = "loglik")

# Rsquared
rsquared(mod)

# get summary
summary(mod)

# get residuals
residuals(mod)

# plotting
plot(mod)

plot(mod, node_attrs = list(
   shape = "rectangle", color = "black",
   fillcolor = "orange", x = 3, y=1:4))

# add correlated error
mod2 <- psem(
  lm(rich ~ cover, data=keeley),
  lm(cover ~ firesev + age, data=keeley),
  lm(firesev ~ age, data=keeley),
  rich %~~% firesev,
  data = keeley)

# get summary
summary(mod2)

# AIC of two models
anova(mod, mod2)

# test mixed models
library(lme4)
library(nlme)

# using lme
shipley_psem <- psem(
  lme(DD ~ lat, random = ~ 1 | site / tree, na.action = na.omit,
  data = shipley),
  lme(Date ~ DD, random = ~ 1 | site / tree, na.action = na.omit,
  data = shipley),
  lme(Growth ~ Date, random = ~ 1 | site / tree, na.action = na.omit,
  data = shipley),
  glmmPQL(Live ~ Growth + (1 | site) + (1 | tree),
  family = binomial(link = "logit"), data = shipley))

summary(shipley_psem)

# using lme4
shipley_psem_lme4 <- psem(
  lmer(DD ~ lat + (1 | site / tree), 
      data = shipley),
  lmer(Date ~ DD + (1 | site / tree), 
      data = shipley),
  lmer(Growth ~ Date + (1 | site / tree),
      data = shipley),
  glmer(Live ~ Growth + (1 | site) + (1 | tree),
        family = binomial(link = "logit"), data = shipley),
  data = shipley)

summary(shipley_psem_lme4)

# multigroup
data(meadows)

meadows$grazed <- factor(meadows$grazed)

meadow_mod <- psem(
  lm(mass ~ elev, data = meadows),
  lm(rich ~ elev + mass, data = meadows),
  data = meadows
)

multigroup(meadow_mod, group = "grazed")