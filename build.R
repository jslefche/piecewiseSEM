#-----BUILD SCRIPT--------------------------------------------------------------

#-----Generate documentation----------------------------------------------------

roxygen2::roxygenise()

#-----Create website for package------------------------------------------------

pkgdown::build_site()

#-----Build package-------------------------------------------------------------

# piecewiseSEM <- devtools::as.package("./piecewiseSEM")

# Add files to .Rbuildignore
usethis::use_build_ignore(c("build.R", ".git", ".gitignore", "docs"))

# Load package
devtools::load_all(".", reset = T)

# Check and build
devtools::check(".", cran = T)

devtools::build(".")

devtools::check_built(".")

# Check on R-hub

# devtools::spell_check()

# devtools::check_rhub("piecewiseSEM")

# check_win_devel()

#-----Examples------------------------------------------------------------------

#run_examples(piecewiseSEM) 

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

AIC(mod) #likelihood based

AIC(mod, AIC.type = "dsep") #dsep based

# Chisq
LLchisq(mod)

AIC(mod, AIC.type = "loglik")

# Rsquared
rsquared(mod)

# get summary
summary(mod)

# get residuals
residuals(mod)

# plotting
# plot(mod)
# 
# plot(mod, node_attrs = list(
#    shape = "rectangle", color = "black",
#    fillcolor = "orange", x = 3, y = 1:4))

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
# library(MASS)
library(nlme)

# using lme
shipley_psem <- psem(
  lme(DD ~ lat, random = ~ 1 | site / tree, na.action = na.omit,
  data = shipley),
  lme(Date ~ DD, random = ~ 1 | site / tree, na.action = na.omit,
  data = shipley),
  lme(Growth ~ Date, random = ~ 1 | site / tree, na.action = na.omit,
  data = shipley),
  glmer(Live ~ Growth + (1 | site) + (1 | tree),
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

#-----Submit to CRAN------------------------------------------------------------

devtools::release("piecewiseSEM")
