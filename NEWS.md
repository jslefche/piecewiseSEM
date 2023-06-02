---
editor_options: 
  markdown: 
    wrap: 72
---

## Version 2.3.01 Release Notes
-   Updated contact email

## Version 2.3.0 Release Notes

-   Now supports model classes: `glmmTMB` , `Sarlm`
-   Added d-sep and Chi-squared tests for `gam` models
-   Added support for `tibble` and `data.table`
-   Added `stop_psem` to deal with unsupported functions, e.g., `poly`
-   Fixed issue with infinite loops with `sortDag`
-   Fixed issue with `getSDx` not being able to get numeric from
    data.frame with additional attributes
-   Fixed issue with `GetData` obtaining transformed data from `lme4`
    objects
-   Fixed issue with non-recursive models and `igraph::is_dag`
-   Fixed issue with `GetSDy` not returning correlation when NAs are
    present
-   Fixed issue with arguments in formulae for `all.vars_trans` and
    `all.vars_notrans`
-   Fixed issue with `GetSDy` and `glm.nb`
-   Fixed issue with `rsquared` when random slopes were not present as
    corresponding fixed effects
-   Fixed issue with `GetSingleData` getting transformed data from `gam`
    models
-   Fixed issue with `anovaLRT` returning weird rownames
-   Fixed issue passing single object to `anova.psem`

## Version 2.2.0 Release Notes

-   Added log-likelihood chi-squared goodness-of-fit and AIC statistics
    from Shipley & Douma 2020
-   Added coefficient standardization for GLM(M)
-   Added preliminary support for GAMs
-   Fixed issue with interactions coded with ":" and `coefs`
-   Fixed issue with AICc and `AIC.type = "dsep"` calculating
    incorrectly

## Version 2.1.2 Release Note

-   Updated `coefs` to remove depracated function `emmeans::CLD`

## Version 2.1.1 Release Notes

-   Fixed bug with interactions without main effects in `coefs`
-   Fixed issue with `rsquared.glm` for "gaussian" family
    and`method = "nagelkerke" -Fixed issue with`car::Anova`output and`handleCategoricalCoefs\`

## Version 2.1.0 Release Notes

-   implemented `igraph` package for DAG construction and evaluation
-   replaced `pbkrtest` with `car::Anova`
-   added `plot.psem` function using `DiagrammR` package
-   `coefs` now returns coefficients for categorical variables using
    `emmeans` and `Anova`
-   added `multigroup` function and analysis
-   added sqrt link for poisson GLMMs
-   Changed contact e-mail to new affiliation
-   Fixed bug with `KRp` and `lmerMod` models with intercepts
-   Fixed bug with `getOLRE` and single observation-level random effects
    glmerMods
-   Fixed errant parentheses in `rsquared` leading to wrong values for
    lme
-   Fixed bug in `cyclic` where incorrect error was returned

## Version 2.0.2 Release Notes

-   Fixed bug with `partialCorr` and negative correlations returning
    wrong P-value
-   Fixed with `coefs` and standardization with mixed models

## Version 2.0.1 Release Notes

-   Fixed bug with `KRp` and uneven sample size
-   New warning issued when NAs present in the dataset
-   Added Gamma distribution to `rsquared`
-   Fixed bug to determine whether graph is cyclic
-   Fixed bug thowing errors with calls to `lmer`
-   Fixed bug creating errors with objects fit using `lmerTest`
-   Added pkgdown website

## Version 2.0 Release Notes

### New syntax

-   All functions have been re-written from the ground up
-   Incorporates new `psem` function and S3 objects
-   All necessary information can now be obtained with a single function
    `summary`

### Updated R[2] functions

-   Extends to new distributions and model types using a single function
    `rsquared`

### New standardization procedures

-   Implements range standardization for all response types
-   Adds multiple forms of standardization for binary responses, see
    `?coefs`

### NOTES

-   Removed `sem.plot` function
-   Removed `sem.lavaan` function
-   `groups=` argument is currently broken but will be fixed in version
    2.1 (see doc)
