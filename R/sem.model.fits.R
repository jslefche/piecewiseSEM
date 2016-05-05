sem.model.fits = function(modelList, aicc = FALSE) {

  # If object is just an individual model, convert to a list
  if(all(class(modelList) != "list")) modelList = list(modelList)

  # Check to see if classes are supported
  if(!all(sapply(modelList, function(i) 
    
    all(class(i) %in% c("lm", "glm", "gls", "pgls", "lme", "lmerMod", "merModLmerTest", "glmerMod")) 
    
  ) ) ) warning("(Pseudo-)R^2s are not yet supported for some model classes!")
  
  # Check if all responses in the model list are the same
  same = length(unique(unlist(sapply(modelList, function(x) all.vars(formula(x))[1])))) == 1
  
  # Apply functions across all models in the model list
  ret = do.call(rbind, lapply(modelList, function(model) {
    
    # Create return data.frame
    ret = data.frame(
      Class = class(model)[1],
      Family = "gaussian",
      Link = "identity",
      N = nobs(model),
      Marginal = NA,
      Conditional = NA
    )
  
    # Get R2 for class == lm
    if(all(class(model) == "lm")) {
      
      # Extract r-squared from model summary
      ret$Marginal = summary(model)$r.squared
      
      # Retrieve AIC(c)
      if(same == TRUE) {
        
        if(aicc == FALSE) 
          
          ret$AIC = AIC(model) else
            
            ret$AICc = AIC(model) + (2 * (attr(logLik(model), "df")) * (attr(logLik(model), "df") + 1)) / (nobs(model) - attr(logLik(model), "df") - 1)
          
      }
        
    }
      
    # Get R2 for class == glm
    if(any(class(model) %in% c("glm", "gls", "pgls"))) {
      
      # Classify model family
      if("glm" %in% class(model)) ret$Family = summary(model)$family[[1]]
      
      # Classify link function
      if("glm" %in% class(model)) ret$Link = summary(model)$family[[2]]
      
      # Calculate R2 values
      if(any(class(model) %in% c("glm"))) ret$Marginal = (1 - (model$deviance / model$null.deviance)) else 
        
        if(any(class(model) %in% c("pgls"))) ret$Marginal = 1 - (model$RSSQ / model$NSSQ) else {
          
          null.mod = update(model, . ~ 1)
          
          ret$Marginal = 1 - (model$sigma / null.mod$sigma)^2
          
          }
  
      # Calculate model AIC
      if(same == TRUE) {
        
        if(aicc == FALSE) 
          
          ret$AIC = AIC(model) else
            
            ret$AICc = AIC(model) + (2 * (attr(logLik(model), "df")) * (attr(logLik(model), "df") + 1)) / (nobs(model) - attr(logLik(model), "df") - 1)
          
      }
    
    }
    
    # Get R2 for class == merMod
    if(any(class(model) %in% c("lmerMod", "merModLmerTest"))) {
      
      # Test for non-zero random effects
      if(any(sapply(VarCorr(model), function(x) !all(attr(x, "stddev") > 0))))
         
         stop("Some variance components equal zero. Respecify random structure!")

      # Get variance of fixed effects by multiplying coefficients by design matrix
      varF = var(as.vector(fixef(model) %*% t(model@pp$X)))
      
      # Check to see if random slopes are present as fixed effects
      ref = ranef(model)
      
      ref.names = ifelse(length(ref) > 1, sapply(ref, names), names(ref))
      
      if(any(!ref.names %in% names(fixef(model))))
        
        stop("Random slopes not present as fixed effects. This artificially inflates calculations of conditional R2. Respecify fixed structure!")
      
      # Separate observation variance from variance of random effects
      n.obs = names(unlist(lapply(ranef(model), nrow))[!unlist(lapply(ranef(model), nrow)) == nrow(model@pp$X)])
      
      # Get variance of random effects 
      varRand = sum(
        
        sapply(VarCorr(model)[n.obs], function(Sigma) { #sapply(VarCorr(model)[n.obs], function(Sigma) {
          
          X = model.matrix(model)  
          
          Z = X[, rownames(Sigma), drop = FALSE]
          
          sum(diag(Z %*% Sigma %*% t(Z)))/nrow(X)
          
        } )
        
      )
      
      # Get residual variance
      varResid = attr(VarCorr(model), "sc")^2
      
      # Calculate R2 values
      ret$Marginal = varF / (varF + varRand + varResid)
        
      ret$Conditional = (varF + varRand) / (varF + varRand + varResid)
      
      # Calculate model AIC
      model.ml = update(model, REML = FALSE)
      
      if(same == TRUE) {
        
        if(aicc == FALSE) 
          
          ret$AIC = AIC(model.ml) else
            
            ret$AICc = AIC(model.ml) + (2 * (attr(logLik(model.ml), "df")) * (attr(logLik(model.ml), "df") + 1)) / (nobs(model.ml) - attr(logLik(model.ml), "df") - 1)
          
      }
      
    }
    
    # Get R2 for class == lme
    if(all(class(model) == "lme")) {
      
      # Test for non-zero random effects
      if(any(sapply(VarCorr(model), function(x) !all(attr(x, "stddev") > 0))))
        
        stop("Some variance components equal zero. Respecify random structure!")
      
      # Get design matrix of fixed effects from model
      Fmat = model.matrix(eval(model$call$fixed)[-2], model$data)
      
      # Remove omitted observations and unused factor levels
      if(!is.null(model$na.action)) Fmat = Fmat[-model$na.action, match(names(fixef(model)), colnames(Fmat))]
      
      # Get variance of fixed effects by multiplying coefficients by design matrix
      varF = var(as.vector(fixef(model) %*% t(Fmat)))

      # Check to see if random slopes are present as fixed effects
      ref = ranef(model)
      
      ref.names = ifelse(length(ref) > 1, sapply(ref, names), names(ref))
      
      if(any(!ref.names %in% names(fixef(model))))
        
        stop("Random slopes not present as fixed effects. This artificially inflates calculations of conditional R2. Respecify fixed structure!")
      
      # Get variance of random effects
      if(any(class(try(getVarCov(model), silent = TRUE)) == "try-error")) {
      
        Sigma = suppressWarnings(as.numeric(VarCorr(model)[, 1]))
      
        names(Sigma) = gsub(" =", "", rownames(VarCorr(model)))
        
        Sigma = as.list(na.omit(Sigma[-length(Sigma)]))
          
        Sigma = lapply(1:length(Sigma), function(i) matrix(Sigma[[i]], dimnames = list(names(Sigma)[i], names(Sigma)[i])) )
        
        varRand = sum(
          
          sapply(Sigma, function(Sigma1) { 
            
            Z = Fmat[, colnames(Sigma1), drop = FALSE]
            
            sum(diag(Z %*% Sigma1 %*% t(Z))) / nrow(Fmat)
            
          } )
          
        )
        
      } else {
        
        Z = Fmat[, rownames(VarCorr(model))[which(rownames(VarCorr(model)) %in% colnames(Fmat))], drop = FALSE]
        
        varRand = sum(diag(Z %*% getVarCov(model) %*% t(Z))) / nrow(Fmat)
        
      }
      
      # Get residual variance
      varResid = summary(model)$sigma^2
      
      # Calculate R2 values
      ret$Marginal = varF / (varF + varRand + varResid)
      
      ret$Conditional = (varF + varRand) / (varF + varRand + varResid)
      
      # Calculate model AIC
      model.ml = update(model, data = model$data, method = "ML")
      
      if(same == TRUE) {
        
        if(aicc == FALSE) 
          
          ret$AIC = AIC(model.ml) else
            
            ret$AICc = AIC(model.ml) + (2 * (attr(logLik(model.ml), "df")) * (attr(logLik(model.ml), "df") + 1)) / (nobs(model.ml) - attr(logLik(model.ml), "df") - 1)
          
      }

    }
    
    # Get R2 for class == "glmerMod"
    if(any(class(model) == "glmerMod")) {
      
      # Test for non-zero random effects
      if(any(sapply(VarCorr(model), function(x) !all(attr(x, "stddev") > 0))))
        
        stop("Some variance components equal zero. Respecify random structure!")
      
      # Classify model family
      ret$Family = summary(model)$family
      
      # Classify link function
      ret$Link = summary(model)$link
      
      # Get variance of fixed effects by multiplying coefficients by design matrix
      varF = var(as.vector(fixef(model) %*% t(model@pp$X)))
      
      # Check to see if random slopes are present as fixed effects
      ref = ranef(model)
      
      ref.names = ifelse(length(ref) > 1, sapply(ref, names), names(ref))
      
      if(any(!ref.names %in% names(fixef(model))))
        
        stop("Random slopes not present as fixed effects. This artificially inflates calculations of conditional R2. Respecify fixed structure!")
      
      # Separate observation variance from variance of random effects
      n.obs = names(unlist(lapply(ranef(model), nrow))[!unlist(lapply(ranef(model), nrow)) == nrow(model@pp$X)])

      # Get variance of random effects 
      varRand = sum(
        
        sapply(VarCorr(model)[n.obs], function(Sigma) {
          
          X = model.matrix(model)  
          
          Z = X[, rownames(Sigma), drop = FALSE]
          
          sum(diag(Z %*% Sigma %*% t(Z)))/nrow(X)
          
        } )
        
      )
      
      # Get overdispersion variance
      obs = names(unlist(lapply(ranef(model), nrow))[unlist(lapply(ranef(model), nrow)) == nrow(model@pp$X)])
      
      if(length(obs) == 0) varDisp = 0 else {
        
        varDisp =  sum(
          
          sapply(VarCorr(model)[obs], function(Sigma) {
            
            X = model.matrix(model)  
            
            Z = X[, rownames(Sigma)]
            
            sum(diag(Z %*% Sigma %*% t(Z)))/nrow(X)
            
          } )
          
        )
        
      }
      
      # Get distribution-specific variance
      if(ret$Family == "binomial") {
        
        if(ret$Link == "logit") varDist = (pi^2)/3
        
        else if(ret$Link == "probit") varDist = 1 else {
          
          warning(paste("Model link '", summary(model)$link, "' is not yet supported for the ", summary(model)$family, "distribution"))
          
          varDist = NA
          
        }
        
      } else if(ret$Family == "poisson") {
        
        # Generate null model (intercept and random effects only, no fixed effects)
        null.model = update(model, formula = paste(". ~ ", get.random.formula(model, "~1", modelList = NULL)))
        
        # Get the fixed effects of the null model
        null.fixef = as.numeric(fixef(null.model))
        
        if(ret$Link == "log") varDist = log(1 + 1/exp(null.fixef))
        
      } else if(ret$Link == "sqrt") varDist = 0.25 else {
        
        warning(paste("Model link '", summary(model)$link, "' is not yet supported for the ", summary(model)$family, "distribution"))
        
        varDist = NA
        
      }
      
      # Calculate R2 values
      ret$Marginal = varF / (varF + varRand + varDisp + varDist)
      
      ret$Conditional = (varF + varRand) / (varF + varRand + varDisp + varDist)
      
      # Calculate model AIC
      if(same == TRUE) {
        
        if(aicc == FALSE) 
          
          ret$AIC = AIC(model) else
            
            ret$AICc = AIC(model) + (2 * (attr(logLik(model), "df")) * (attr(logLik(model), "df") + 1)) / (nobs(model) - attr(logLik(model), "df") - 1)
          
      }
      
    }
    
#     # Get R2 for class == "glmmPQL"
#     if(any(class(model) == "glmmPQL")) {
#       
#       # Classify model family
#       ret$Family = summary(model)$family[1]
#       
#       # Classify link function
#       ret$Link = summary(model)$family[2]
#       
#       # Get design matrix of fixed effects from model
#       Fmat = model.matrix(eval(model$call$fixed)[-2], model$data)
#       
#       # Remove omitted observations and unused factor levels
#       if(!is.null(model$na.action)) Fmat = Fmat[-model$na.action, match(names(fixef(model)), colnames(Fmat))]
#    
#       # Get variance of fixed effects by multiplying coefficients by design matrix
#       varF = var(as.vector(fixef(model) %*% t(Fmat)))
#       
#       # Get random effects design matrix
#       rand.mat = ranef(model)
#       
#       # If not in a list, convert and name levels
#       if(!any(class(rand.mat) %in% "list")) {
#         
#         rand.mat = list(rand.mat)
#         
#         names(rand.mat) = attr(unclass(getVarCov(model)), "group.levels")
#         
#       }
#       
#       # Separate observation variance from variance of random effects
#       n.obs = names(rand.mat)[!sapply(rand.mat, nrow) == nrow(Fmat)]
#       
#       Sigma = if(length(n.obs) == 1) as.numeric(VarCorr(model)[1]) else
#         
#         as.numeric(VarCorr(model)[match(paste(n.obs, "="), rownames(VarCorr(model))) + 1, 1])
#       
#       names(Sigma) = if(length(n.obs) == 1) rownames(VarCorr(model))[1] else 
#         
#         rownames(VarCorr(model))[which(rownames(VarCorr(model)) %in% paste(n.obs, "=")) + 1]
# 
#       # Get variance of random effects, excluding observation variance
#       varRand = sum(
#         
#         sapply(Sigma, function(Sigma1) { 
#           
#           X = model.matrix(model)  
#           
#           Z = Fmat[, names(Sigma1), drop = FALSE]
#           
#           sum(diag(Z %*% Sigma1 %*% t(Z))) / nrow(Fmat)
#           
#         } )
#         
#       )
#     
#       # Get overdispersion variance
#       obs = names(rand.mat)[sapply(rand.mat, nrow) == nrow(Fmat)]
# 
#       if(length(obs) == 0) varDisp = 0 else {
#         
#         varDisp =  sum(
#           
#           sapply(VarCorr(model)[obs], function(Sigma) {
#             
#             Z = Fmat[, rownames(Sigma), drop = FALSE]
#             
#             sum(diag(Z %*% Sigma %*% t(Z))) / nrow(Fmat)
#             
#           } )
#           
#         )
#         
#       }
#       
#       # Get distribution-specific variance
#       if(ret$Family == "binomial") {
#         
#         if(ret$Link == "logit") varDist = (pi^2)/3
#         
#         else if(ret$Link == "probit") varDist = 1 else {
# 
#           warning(paste("Model link '", summary(model)$family[2], "' is not yet supported for the ", summary(model)$family[1], "distribution"))
#           
#           varDist = NA
#           
#         }
#         
#       } else if(ret$Family == "poisson") {
#         
#         # Generate null model (intercept and random effects only, no fixed effects)
#         null.model = update(model, fixed = ". ~ 1")
#         
#         # Get the fixed effects of the null model
#         null.fixef = as.numeric(fixef(null.model))
#         
#         if(ret$Link == "log") varDist = log(1 + 1/exp(null.fixef))
#         
#       } else if(ret$Link == "sqrt") varDist = 0.25 else {
#         
#         warning(paste("Model link '", summary(model)$family[2], "' is not yet supported for the ", summary(model)$family[1], "distribution"))
#         
#         varDist = NA
#         
#       }
#       
#       # Calculate R2 values
#       ret$Marginal = varF / (varF + varRand + varDisp + varDist)
#       
#       ret$Conditional = (varF + varRand) / (varF + varRand + varDisp + varDist)
#       
#       # Calculate model AIC
#       ret$AIC = NA
#       
#     }
    
    # Return results
    return(ret)
    
  } ) )
  
  if(any(ret$N <= 40) & colnames(ret)[ncol(ret)] != "AICc") warning("N < 40, consider using aicc = TRUE")
  
  # Get list of response vectors
  resp = sapply(modelList, function(x) all.vars(formula(x))[1])
  
  # Calculate delta AIC
  if(length(resp) > 1 & all(resp[1] == resp) & any(!is.na(ret[, ncol(ret)]))) {
     
    ret = cbind(
      ret, 
      ret[, ncol(ret)] - min(ret[, ncol(ret)])
    )
    
    colnames(ret)[ncol(ret)] = paste0("d", colnames(ret)[ncol(ret) - 1]) #intToUtf8(0x0394), colnames(ret)[ncol(ret) - 1])
    
  }
  
  return(ret)
  
}