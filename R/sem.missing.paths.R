sem.missing.paths = function(
 
  modelList, data, corr.errors = NULL, add.vars = NULL, grouping.vars = NULL, top.level.vars = NULL, filter.exog = TRUE,
  adjust.p = FALSE, basis.set = NULL, model.control = NULL, .progressBar = TRUE
  
  ) {
  
  # Get basis set
  if(is.null(basis.set)) basis.set = sem.basis.set(modelList, corr.errors, add.vars)   

  # Filter exogenous variables
  if(filter.exog == TRUE) basis.set = filter.exogenous(modelList, basis.set, corr.errors, add.vars) 

  # Add progress bar
  if(.progressBar == T) pb = txtProgressBar(min = 0, max = length(basis.set), style = 3) else pb = NULL
  
  # Perform d-sep tests
  pvalues.df = do.call(rbind, lapply(1:length(basis.set), function(i) {
    
    # Get basis model from which to build the d-sep test
    basis.mod = modelList[[match(basis.set[[i]][2], sapply(modelList, function(j) strsplit(paste(formula(j)), "~")[[2]]))]]
    
    # Get fixed formula
    rhs = if(length(basis.set[[i]]) <= 2) paste(basis.set[[i]][1]) else
      
      paste(basis.set[[i]][c(1,3:length(basis.set[[i]]))], collapse = " + ")
    
    # Get random formula
    random.formula = get.random.formula(basis.mod, rhs, modelList)
    
    # Aggregate at the level of the grouping variable
    if(!is.null(grouping.vars) & any(top.level.vars %in% gsub(".*\\((.*)\\)", "\\1", unlist(as.character(formula(basis.mod)[2])))))
      
      data = suppressWarnings(aggregate(data, by = lapply(grouping.vars, function(i) data[ ,i]), mean, na.rm = T))
    
    # Get model controls
    control = get.model.control(basis.mod, model.control)
    
    # Update basis model with new formula and random structure based on d-sep
    basis.mod = suppressWarnings(if(is.null(random.formula)) 
      
      update(basis.mod, formula(paste(basis.set[[i]][2], " ~ ", rhs)), data = data) else
        
        if(any(class(basis.mod) %in% c("lme", "glmmPQL"))) 
          
          update(basis.mod, fixed = formula(paste(basis.set[[i]][2], " ~ ", rhs)), random = formula(random.formula), control = control, data = data) else
            
            update(basis.mod, formula = formula(paste(basis.set[[i]][2], " ~ ", rhs, " + ", random.formula, sep = "")), control = control, data = data) 
    )
    
    if(any(class(basis.mod) %in% "lmerMod")) basis.mod = as(basis.mod, "merModLmerTest") 
    
    # Insert stop if lmerTest does not converge
    if(all(class(basis.mod) %in% c("lmerMod", "merModLmerTest"))) {
      
      x = try(suppressMessages(suppressWarnings(summary(basis.mod))$coefficients[1,5]), silent = T)
      
      if(class(x) == "try-error") stop("lmerTest did not converge, no p-values to report. Consider specifying model.control or using nlme package") 
      
    }
    
    # Get row number from coefficient table for d-sep variable
    row.num = which(basis.set[[i]][1] == attr(terms(basis.mod), "term.labels")) + 1
    
    # Get row number if interaction variables are switched
    if(length(row.num) == 0 & grepl("\\:", basis.set[[i]][1]))
      
      row.num = which(paste(rev(strsplit(basis.set[[i]][1], ":")[[1]]), collapse = ":") == attr(terms(basis.mod), "term.labels")) + 1
    
    # Return new coefficient table
    ret = if(any(class(basis.mod) %in% c("lm", "glm", "negbin", "pgls", "glmerMod", "merModLmerTest")))
      
      as.data.frame(t(unname(summary(basis.mod)$coefficients[row.num, ]))) else
      
        as.data.frame(t(unname(summary(basis.mod)$tTable[row.num, ])))  
    
    # Add df if summary table does not return
    if(length(ret) != 5 & any(class(basis.mod) %in% c("lm", "glm", "negbin", "pgls"))) 
      
      ret = cbind(ret[1:2], summary(basis.mod)$df[2], ret[3:4]) else
        
        if(length(ret) != 5)
          
          ret = cbind(ret[1:2], "", ret[3:4])
      
    # Rename columns 
    names(ret) = c("estimate", "std.error", "DF", "crit.value", "p.value")
    
    # Adjust p-value based on Shipley 2013
    if(adjust.p == TRUE) {
      
      if(all(class(basis.mod) %in% c("lme", "glmmPQL"))) {
        
        t.value = summary(basis.mod)$tTable[row.num, 4] 
        
        ret[5] = 2*(1 - pt(abs(t.value), nobs(basis.mod) - sum(apply(basis.mod$groups, 2, function(x) length(unique(x))))))
        
      } else if(all(class(basis.mod) %in% c("glmerMod", "merModLmerTest"))) {
        
        z.value = summary(basis.mod)$coefficients[row.num, 4]
        
        ret[5] = 2*(1 - pt(abs(z.value), nobs(basis.mod) - sum(summary(basis.mod)$ngrps))) } 
      
    }
    
    if(.progressBar == TRUE) setTxtProgressBar(pb, i)
    
    # Bind in d-sep metadata
    cbind(missing.path = paste(basis.set[[i]][2], "<-", paste(basis.set[[i]][1], collapse = "+")), ret)
    
  } ) )
  
  if(!is.null(pb)) close(pb)  
  
  return(pvalues.df)
  
}