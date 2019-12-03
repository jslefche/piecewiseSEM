###### Build script
library(devtools)
setwd("../")
# library(Rd2roxygen)
# Rd2roxygen::Rd2roxygen("piecewiseSEM")


piecewiseSEM <- as.package("./piecewiseSEM") 
devtools::use_build_ignore(c("build.R", ".git", ".gitignore", "docs"),
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
#run_examples(piecewiseSEM) 

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
basisSet(mod)
dSep(mod)

summary(mod)

res <- residuals(mod)
anova(mod)
fisherC(mod)

plot(mod)

plot(mod, node_attrs = list(
   shape = "rectangle", color = "black",
   fillcolor = "orange", x = 3, y=1:4))

mod2 <- psem(
  lm(rich ~ cover, data=keeley),
  lm(cover ~ firesev + age, data=keeley),
  lm(firesev ~ age, data=keeley),
  rich %~~% firesev,
  data = keeley
  
)
summary(mod2)

anova(mod, mod2)

#######################################

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

#testing interactions
data(meadows)
meadows$grazed <- factor(meadows$grazed)

meadow_mod <- psem(
  lm(mass ~ grazed*elev, data = meadows),
  lm(rich ~ grazed*(elev + mass), data = meadows),
  data = meadows
)

dSep(meadow_mod)

anova(meadow_mod)

meadom_mod_const <- psem(
  lm(mass ~ grazed + elev, data = meadows),
  lm(rich ~ grazed +  mass, data = meadows),
  data = meadows
)


anova(meadom_mod_const)


anova(meadow_mod, meadom_mod_const)

fisherC(meadow_mod)
fisherC(meadom_mod_const)


#make sure factors are included in conditional independence tests
meadows$grazed <- factor(meadows$grazed)
meadom_mod_graz <- psem(
  lm(mass ~  elev, data = meadows),
  lm(rich ~ grazed +  mass, data = meadows),
  data = meadows
)

dSep(meadom_mod_graz)

#

#lmerTest
# library(lmerTest)
# shipley_psem_lmerTest <- psem(
#   
#   lmer(DD ~ lat + (1 | site / tree), 
#        data = shipley),
#   
#   lmer(Date ~ DD + (1 | site / tree), 
#        data = shipley),
#   
#   lmer(Growth ~ Date + (1 | site / tree),
#        data = shipley),
#   
#   glmer(Live ~ Growth + (1 | site) + (1 | tree),
#         family = binomial(link = "logit"), data = shipley),
#   
#   data = shipley
#   
# )
# 
# summary(shipley_psem_lmerTest)
