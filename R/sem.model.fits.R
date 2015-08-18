sem.model.fits = function(modelList) {

  # If object is just an individual model, convert to a list
  if(all(class(modelList) != "list")) modelList = list(modelList)

  # Check to see if classes are supported
  if(!all(sapply(modelList, function(i) 
    
    all(class(i) %in% c("lm", "glm", "negbin", "gls", "pgls", "lme", "lmerMod", "merModLmerTest", "glmerMod", "glmmPQL")) 
    
  ) ) ) warning("Pseudo-R2s are not yet supported for some model classes!")
  
  # Apply functions across all models in the model list
  do.call(rbind, lapply(modelList, function(model) {
    
    # Create return data.frame
    ret = data.frame(
      Class = class(model)[1],
      Family = "gaussian",
      Link = "identity",
      Marginal = NA,
      Conditional = NA,
      AIC = AIC(model) 
    )
  
    # Get R2 for class == lm
    if(all(class(model) == "lm")) ret$Marginal = summary(model)$r.squared
      
    # Get R2 for class == glm
    if(any(class(model) %in% c("glm", "gls", "pgls"))) {
      
      # Classify model family
      if("glm" %in% class(model)) ret$Family = summary(model)$family[[1]]
      
      # Classify link function
      if("glm" %in% class(model)) ret$Link = summary(model)$family[[2]]
      
      # Calculate McFadden R2 values
      ret$Marginal = 1 - (as.numeric(logLik(model)) / as.numeric(logLik(update(model, ". ~ 1"))))

      ret$Conditional = NA
                    
      # Calculate model AIC
      ret$AIC = AIC(model)
      
    }
    
    # Get R2 for class == merMod
    if(any(class(model) %in% c("lmerMod", "merModLmerTest"))) {
      
      # Get variance of fixed effects by multiplying coefficients by design matrix
      varF = var(as.vector(fixef(model) %*% t(model@pp$X)))
   
      # Separate observation variance from variance of random effects
      n.obs = names(unlist(lapply(ranef(model), nrow))[!unlist(lapply(ranef(model), nrow)) == nrow(Fmat)])
      
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
      ret$AIC = AIC(update(model, REML = FALSE))
      
    }
    
    # Get R2 for class == lme
    if(all(class(model) == "lme")) {
      
      # Get design matrix of fixed effects from model
      Fmat = model.matrix(eval(model$call$fixed)[-2], model$data)
      
      # Get variance of fixed effects by multiplying coefficients by design matrix
      varF = var(as.vector(fixef(model) %*% t(Fmat)))
 
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
      ret$AIC = AIC(update(model, method = "ML"))
      
    }
    
    # Get R2 for class == "glmerMod"
    if(any(class(model) == "glmerMod")) {
      
      # Classify model family
      ret$Family = summary(model)$family
      
      # Classify link function
      ret$Link = summary(model)$link
      
      # Get variance of fixed effects by multiplying coefficients by design matrix
      varF = var(as.vector(fixef(model) %*% t(model@pp$X)))
      
      # Separate observation variance from variance of random effects
      n.obs = names(unlist(lapply(ranef(model), nrow))[!unlist(lapply(ranef(model), nrow)) == nrow(Fmat)])
      
      # Get variance of random effects 
      varRand = sum(
        
        sapply(VarCorr(model)[n.obs], function(Sigma) {
          
          X = model.matrix(model)  
          
          Z = X[, rownames(Sigma), drop = FALSE]
          
          sum(diag(Z %*% Sigma %*% t(Z)))/nrow(X)
          
        } )
        
      )
      
      # Get overdispersion variance
      obs = names(unlist(lapply(ranef(model), nrow))[unlist(lapply(ranef(model), nrow)) == nrow(Fmat)])
      
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
      ret$AIC = AIC(model)
      
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
      
}
