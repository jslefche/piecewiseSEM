# piecewiseSEM Change Log

### 2016-09-19 Version 1.2
  * Feature addition: new function `acyclic` to test for acyclic DAGs
  * Feature addition: new function `rsquared` for easier calling of `sem.model.fits`
  * Bug fix: order of adjacency matrix in `get.sort.dag`
  * Feature addition: model class `glmmTMB`
  * Bug fix: inclusion of transformed variables as additional boxes in `sem.plot`
  * Bug fix: specification of invalid models in `endogenous.reverse`
  * Feature addition: implemented two-column binding as response in GLM(M)s
  * Bug fix: memory allocation error in `sem.model.fits`
  * Bug fix: inclusion of entries in the basis set without corresponding models in the model list in `endogenous.reverse`
  * Bug fix: crossed random effects in `get.random.formula` not returning correct random slopes
  * Bug fix: `gls` models in `sem.missing.paths` could not update basis model

### 2016-07-05 Version 1.1.3
  * Feature addition: report of intercepts in `sem.coefs`
  * Bug fix: warning with random slopes not present as fixed effects in `sem.model.fits`
  * Feature addition: new function `endogenous.reverse` to assist with independence claims among intermediate endogenous variables not fit to a normal distribution
  * Bug fix: negative binomial models fail to return R2s in `sem.model.fits`
  * Bug fix: Failed to index correct P-value in `sem.missing.paths` when categorical variables are present
  * Bug fix: intermediate endogenous variables not properly removed from the basis set in `sem.missing.paths`

### 2016-06-06 Version 1.1.2
  * Feature addition: significance indicators for P-value outputs
  * Bug fix: remove random effects from `get.scaled.data`
  * Feature addition: Use `KRmodcomp` from `pbkrtest` instead of `get_LB_dff` to calculate P-values for `lmer` models
  * Bug fix: transformations for scaled coefficients no longer fails in `sem.coefs`
  * Bug fix: all interactions now show up in `sem.coefs` when scaled
  * Bug fix: random slopes as fixed effects not registering properly in `sem.model.fits`
  * Bug fix: added \dontrun{} to all help files

### 2016-05-18 Version 1.1.1
  * Reduced time it takes examples to run
  * Stop: if duplicate responses are detected in the model list
  * Temporary bug fix: returns lowest P-value for when d-sep tests include non-linear intermediate endogenous variables, only if family is not Gaussian
  * Bug fix: reversed order of conditioning variables so they come first 
  * Bug fix: removed entries from the basis set that attempt to predict an interaction that does not appear in the model list

### 2016-03-10 Version 1.1.0
  * Feature addition: rudimentary plotting using `sem.plot`
  * Feature addition: support for `glmmadmb` models
  * Feature addition: `pbkrtest` to reliably extract p-values from `lme4` models
  * *Major Bug Fix:* incorrect basis set specified by `ggm::topSort`; introduced `sort.dag` function, removed dependency on `ggm`
  * *Major Bug Fix:* `sem.basis.set` now treats transformed variables as untransformed, resolving duplicate vars in the basis set
  * Bug fix: resolved switched interactions (e.g., x1:x2 vs x2:x1) leading to duplicate vars in the basis set using `get.dag`
  * Bug fix: fixed long standing bug with `lmer` models not returning p-values
  
### 2016-02-09 Version 1.0.4
  * Bug fix: correlated errors among exogenous variables in `sem.coefs`
  * Bug fix: correctly scale data from `pgls` models using `get.scaled.data`
  * Bug fix: corrected basis set for models created with `gls`
  * Feature addition: helper function `get.scaled.model` to get model from scaled coefficients
  * Feature addition: (partial) correlations in sem.coefs
  * Bug fix: passing of correlated errors to lavaan
  * Bug fix: standardized coefficients for variables transformed in model formula

### 2016-01-20 Version 1.0.3
  * Bug fix: get.model.control updated for latest versions of lme4 and nlme
  * Bug fix: hand compute interactions in `sem.lavaan` using argument `compute.int = TRUE`
  * Feature addition: function `get.scaled.data` to handle transformed variables in `sem.coefs` when standardize = "scale"

### 2016-01-15 Version 1.0.2
  * Feature addition: added additional plotting arguments for partial.residuals
  * Feature addition: AICc and delta AIC to sem.model.fits
  * Bug fix: Offsets in model formula treated as predictors in d-sep tests

### 2015-12-13 Version 1.0.1
  * Feature addition: AICc and delta AIC to sem.model.fits
  * Bug fix: issue with interactions in basis set and d-sep tests in sem.missing.paths
  * Bug fix: issue with fixed intercept models and calculating the basis set
  * Bug fix: issue with offset() variables in sem.missing.paths
  * Bug fix: transformed responses for partial residuals
  * Bug fix: duplicated values in basis set (function DAG in ggm package)

### 2015-10-26 Version 1.0.0
  * First release to CRAN

### 2015-10-23 Version 0.9.9
  * Added vignette
  * Fixed issue with design matrix including omitted observations and dropped levels in sem.model.fits
  * Incorporated interactions into partial.resid
  * Fixed bug with lme residuals in partial.resid and predict.sem
  * Modified sem.coefs to return NA for standardized interactions

### 2015-09-08 Version 0.9.8
  * Added functions get.dag and get.basis.set
  * Fixed issue with interactions in sem.missing.paths
  * Fixed issue in sem.missing.paths where only 1 missing path returned NA instead of p-value
  * Removed argument filter.exog = TRUE
  * Fixed error for add.vars in get.dag
  * Fixed issue with transformed variables and sem.coefs
  * Fixed issue with lme4 models and sem.model.fits
  * Fixed issue with incorrect independence claims in get.basis.set
  * Fixed (nlme) merging error in partial.resid
  * Fixed issue with glmmPQL and get.random.structure
  * Added new citation (arXiv)

### 2015-07-31 Version 0.9.1
  * Implemented workaround in sem.basis.set where output changed based on order of variables in formula with interactions
  * Saturated models now return AIC (and AICc) values (with warning)
  * Added standard errors on predictions for mixed models in predict.sem()
  * Fixed typo in output for sem.fisher.c with incorrect df
  * Fixed error with lme4 and gls models in partial.resid()
  * Fixed error with lmerTest returning "merMod" objects
  * New function predict.sem returns model predictions
  * Fixed "unsupported model class" error for merModTest

### 2015-06-15 Version 0.9
  * Major revisions and annotations to all functions to improve transparency and efficiency
  * Included new helper functions get.random.formula() and get.model.control()
  * Now reports df for all model types except glmer
  * Included support for "gls" and phylogenetically independent contrast ("pgls") models
  * Fixed deparse error in get.random.struture()
  * Improved handling of gls() models

### 2015-04-27 Version 0.4.4
  * Fixed issue with transformed corr.errors in get.basis.set
  * get.sem.coefs now does range standardization
  * get.sem.coefs now returns single table with significance values
  * Improved handling of uneven observations in get.partial.resid
  * get.partial.resid now accepts a single model in addition to a list of models
  * get.lavaan.sem now supports additional arguments from lavaan

### 2015-01-26 Version 0.4.3
  * Added new function to return R^2 and AIC values for component models in SEM
  * get.partial.resid now returns residuals plot with fitted line

### 2015-01-03 Version 0.4.2
  * get.missing.paths now returns more information for d-sep tests
  * get.basis.set drops independence claims in the basis set where interactions predict their corresponding main effects 
  * get.fisher.c now reports k df and get.aic now reports K and n df
  * Fixed rounding error that misreported Fisher's C
  * Allow user-specified rounding using `sig=` argument in get.sem.fit, get.fisher.c, and get.aic

### 2014-12-19 Version 0.4.1
  * Fixed minor bug in get.partial.resid where random effects threw an error for lme models
  * Fixed interaction bug in filter.exogenous which was returning incompatible models from the basis set
  * Fixed bug in get.sem.coefs where missing values returned NA for corr.errors
  * Fixed issue with transformed/untransformed variables predicting one another in the basis set
  * Fixed interaction bug in get.missing.paths, filter.exogenous, get.partial.resiod
  * Now can supply single model control or list

### 2014-11-12 Version 0.4.0
  * Fixed bug for transformed variables in get.partial.resid
  * Can now specify an equation as input into get.partial.resid
  * Fixed bug for transformed interaction terms in structured equations for get.basis.set
  * Fixed bug for transformed variables when also specifying raw variables in corr.error
  * Added model control arguments to get.sem.fit and get.partial.resid

### 2014-09-22 Version 0.3.2
  * Added optional switch to display conditional variables in get.missing.paths
  * Removed logLik df from get.aic output

### 2014-09-19 Version 0.3.1
  * Incorporated multilevel data into get.sem.fit
  * Fixed error when length(basis.set) = 1 and .progressBar = T
  * Expanded interaction term in get.basis.set to include 20 letters instead of 10, and run them backwards
  * Renamed variable names in filter.exogenous to better reflect what they are

### 2014-09-12 Version 0.2.6
  * Incorporated transformed variables and multiple interactions into get.sem.coefs
  * Fixed but with add.vars; many functions now require the user to supply a data.frame

### 2014-09-10 Version 0.2.5
  * Standardized argument order across functions
  * Improved handling of correlated errors in get.sem.coefs
  * Incorporated correlated errors into get.sem.lavaan

### 2014-09-05 Version 0.2.3
  * Fixed return of standardized coefficients in get.sem.coefs

### 2014-08-25 Version 0.2.2
  * Added argument to define and quantify correlated errors 
  * Added new function get.partial.resid

### 2014-08-18 Version 0.1.5
  * Added new function get.basis.set (replaced dag.updated)
  * Improved handling of interactions in the basis set
  * Improved handling of p-value rownames in get.missing.paths
  * get.aic now returns df

### 2014-08-04 Version 0.1.4
  * Added new functions get.aic, get.fisher.c, get.sem.coefs, get.lavaan.sem
  * Added README.md
  * Completed all existing .Rd files
  * Updated .progressBar as to not conflict with existing function progressBar
  * Altered handling of model formulae to collapse to a single character string

### 2014-08-03 Version 0.1.0
  * Initial build of package
  * Added CHANGELOG