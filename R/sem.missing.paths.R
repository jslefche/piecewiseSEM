sem.missing.paths = function(
 
  modelList, data, conditional = FALSE, corr.errors = NULL, add.vars = NULL, grouping.vars = NULL, 
  top.level.vars = NULL, adjust.p = FALSE, basis.set = NULL, model.control = NULL, .progressBar = TRUE
  
  ) {
  
  # Get basis set
  if(is.null(basis.set)) basis.set = suppressWarnings(sem.basis.set(modelList, corr.errors, add.vars))

  # Filter exogenous variables
  basis.set = filter.exogenous(modelList, basis.set, corr.errors, add.vars) 

  # Add progress bar
  if(.progressBar == T & length(basis.set) > 0) 
    
    pb = txtProgressBar(min = 0, max = length(basis.set), style = 3) else pb = NULL
  
  # Perform d-sep tests
  if(length(basis.set) > 0) pvalues.df = do.call(rbind, lapply(1:length(basis.set), function(i) {
    
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
    basis.mod.new = suppressWarnings(if(is.null(random.formula)) 
      
      update(basis.mod, formula(paste(basis.set[[i]][2], " ~ ", rhs)), data = data) else
        
        if(any(class(basis.mod) %in% c("lme", "glmmPQL"))) 
          
          update(basis.mod, fixed = formula(paste(basis.set[[i]][2], " ~ ", rhs)), random = random.formula, control = control, data = data) else
            
            update(basis.mod, formula = formula(paste(basis.set[[i]][2], " ~ ", rhs, " + ", random.formula, sep = "")), control = control, data = data) 
      
    )
    
    # Stop if lmerTest does not return p-values
    if("lmerMod" %in% class(basis.mod.new)) {
      
      basis.mod.new = as(basis.mod.new, "merModLmerTest") 
      
      x = try(suppressMessages(summary(basis.mod.new)$coefficients[1,5]), silent = T)
      
      if(class(x) == "try-error") stop("lmerTest did not return p-values. Consider using functions in the nlme package.")
      
    }
    
    # Get row number from coefficient table for d-sep variable
    if(any(!class(basis.mod.new) %in% "pgls")) {
      
      row.num = which(basis.set[[i]][1] == rownames(attr(terms(basis.mod.new), "factors"))[-1]) + 1 
      
      # Get row number if interaction variables are switched
      if(length(row.num) == 0 & grepl("\\:", basis.set[[i]][1])) {
        
        # Get all combinations of interactions
        ints = strsplit(basis.set[[i]][1], ":")
        
        all.ints = sapply(ints, function(x) { 
          
          datf = expand.grid(rep(list(x), length(x)), stringsAsFactors = FALSE)
          
          datf = datf[apply(datf, 1, function(x) !any(duplicated(x))), ]
          
          apply(datf, 1, function(x) paste(x, collapse = ":"))
          
        } )
        
        row.num = which(attr(terms(basis.mod.new), "term.labels") %in% all.ints) + 1
          
        }
      
      } else {
        
        row.num = which(basis.set[[i]][1] == basis.mod.new$varNames)
         
      }
    
    # Return new coefficient table
    ret = if(any(class(basis.mod.new) %in% c("lm", "glm", "negbin", "pgls", "glmerMod", "merModLmerTest")))
      
      as.data.frame(t(unname(summary(basis.mod.new)$coefficients[row.num, ]))) else
      
        as.data.frame(t(unname(summary(basis.mod.new)$tTable[row.num, ])))  
    
    # Add df if summary table does not return
    if(length(ret) != 5 & any(class(basis.mod.new) %in% c("lm", "glm", "negbin", "pgls"))) 
      
      ret = cbind(ret[1:2], summary(basis.mod.new)$df[2], ret[3:4]) else
        
        if(length(ret) != 5)
          
          ret = cbind(ret[1:2], NA, ret[3:4])
      
    # Rename columns 
    names(ret) = c("estimate", "std.error", "DF", "crit.value", "p.value")
    
    # Adjust p-value based on Shipley 2013
    if(adjust.p == TRUE) {
      
      if(all(class(basis.mod.new) %in% c("lme", "glmmPQL"))) {
        
        t.value = summary(basis.mod.new)$tTable[row.num, 4] 
        
        ret[5] = 2*(1 - pt(abs(t.value), nobs(basis.mod.new) - sum(apply(basis.mod.new$groups, 2, function(x) length(unique(x))))))
        
      } else if(all(class(basis.mod.new) %in% c("glmerMod", "merModLmerTest"))) {
        
        z.value = summary(basis.mod.new)$coefficients[row.num, 4]
        
        ret[5] = 2*(1 - pt(abs(z.value), nobs(basis.mod.new) - sum(summary(basis.mod.new)$ngrps))) } 
      
    }
  
    if(.progressBar == TRUE) setTxtProgressBar(pb, i)
    
    # Modify rhs if number of characters exceeds 20
    if(conditional == FALSE) if(nchar(rhs) > 30) rhs = paste(gsub(".\\+.*$", "", rhs), "+ ...")
    
    # Bind in d-sep metadata
    data.frame(missing.path = paste(basis.set[[i]][2], " ~ ", rhs, sep = ""), ret)
    
  } ) ) else
    
    pvalues.df = data.frame(missing.path = NA, estimate = NA, std.error = NA, DF = NA, crit.value = NA, p.value = NA)
  
  # Set degrees of freedom as numeric
  pvalues.df$DF = as.numeric(pvalues.df$DF)
  
  # Output warning if missing.paths contains ...
  if(any(grepl("\\.\\.\\.", pvalues.df$missing.path))) print("Conditional variables have been omitted from output table for clarity")
  
  if(!is.null(pb)) close(pb)  
  
  return(pvalues.df)
  
}