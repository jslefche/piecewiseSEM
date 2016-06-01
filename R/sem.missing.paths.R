sem.missing.paths = function(
 
  modelList, data, conditional = FALSE, corr.errors = NULL, add.vars = NULL, grouping.vars = NULL, 
  adjust.p = FALSE, basis.set = NULL, model.control = NULL, .progressBar = TRUE
  
  ) {
  
  # Get basis set
  if(is.null(basis.set)) basis.set = suppressWarnings(sem.basis.set(modelList, corr.errors, add.vars))

  # Add progress bar
  if(.progressBar == T & length(basis.set) > 0) pb = txtProgressBar(min = 0, max = length(basis.set), style = 3) else pb = NULL
  
  # Identify d-sep tests among endogenous variables and reverse if family != "gaussian"
  names(basis.set) = 1:length(basis.set)
  
  # Get sorted adjacency matrix
  amat = get.sort.dag(get.formula.list(modelList))
  
  # Identify intermediate endogenous variables
  idx = colnames(amat)[amat[colnames(amat)[colSums(amat) == 0], ] > 0]
  
  # Identify variables in the basis set where intermediate endogenous variables are the response
  if(any(sapply(modelList[sapply(modelList, function(x) all.vars(formula(x))[1] %in% idx)], function(x) 
    
    if(any(class(x) %in% c("glm", "negbin", "glmmPQL", "glmerMod"))) FALSE else class(x) != "gaussian"
    
    ))) {
    
    # Add flag
    rev = TRUE
  
    basis.set = append(basis.set, lapply(which(sapply(basis.set, function(i) i[2] %in% idx)), function(i) 
      
      c(basis.set[[i]][2], basis.set[[i]][1], basis.set[[i]][-(1:2)])
      
    ) )
    
  }
  
  # Perform d-sep tests
  if(length(basis.set) > 0) pvalues.df = do.call(rbind, lapply(1:length(basis.set), function(i) {
    
    # Get basis model from which to build the d-sep test
    basis.mod = modelList[[which(sapply(modelList, function(j) {
      
      if(any(class(j) == "pgls")) j = j$formula
      
      rownames(attr(terms(j), "factors"))[1] == basis.set[[i]][2]
    
    } ) ) ]]
      
    # Get fixed formula
    rhs = if(length(basis.set[[i]]) <= 2) paste(basis.set[[i]][1]) else
      
      paste(basis.set[[i]][c(3:length(basis.set[[i]]), 1)], collapse = " + ")
    
    # Get random formula
    random.formula = get.random.formula(basis.mod, rhs, modelList)
    
    # Aggregate at the level of the grouping variable
    if(!is.null(grouping.vars)) {
      
      # Test to see if response is identical for levels of grouping factor
      response = as.character(formula(basis.mod)[2])
      
      response.test = by(data, lapply(grouping.vars, function(j) data[ ,j]), function(x) length(unique(x[, response])) == 1)
       
      response.test = na.omit(vapply(response.test, unlist, unlist(response.test[[1]])))
       
      # If so, aggregate by grouping.vars
      if(all(response.test == TRUE)) {
        
        # Get named list
        groups = lapply(grouping.vars, function(j) data[ ,j])
        
        names(groups) = grouping.vars
        
        # Aggregate and replace data
        data = suppressWarnings(
          
          aggregate(data, by = groups, mean, na.rm = T)
          
          )
        
        # Remove duplicated colnames (keep first instance, from above)
        data = data[, !duplicated(colnames(data))]
        
      }
    
    }
    
    # Get model controls
    control = get.model.control(basis.mod, model.control)
    
    # Update basis model with new formula and random structure based on d-sep
    basis.mod.new = suppressWarnings(
      
      if(is.null(random.formula) | class(basis.mod) == "glmmadmb") 
      
        update(basis.mod, formula(paste(basis.set[[i]][2], " ~ ", rhs)), data = data) else
        
          if(any(class(basis.mod) %in% c("lme", "glmmPQL"))) 
          
            update(basis.mod, fixed = formula(paste(basis.set[[i]][2], " ~ ", rhs)), random = random.formula, control = control, data = data) else
            
              update(basis.mod, formula = formula(paste(basis.set[[i]][2], " ~ ", rhs, " + ", random.formula, sep = "")), control = control, data = data) 
      
    )

    # Get row number from coefficient table for d-sep variable
    # if(any(!class(basis.mod.new) %in% c("pgls"))) {
    #   
    #   # Get row number of d-sep claim
    #   row.num = which(basis.set[[i]][1] == rownames(attr(terms(basis.mod.new), "factors"))[-1]) + 1 
    #   
    #   # Get row number if interaction variables are switched
    #   if(length(row.num) == 0 & grepl("\\:|\\*", basis.set[[i]][1])) {
    #     
    #     # If interaction is reported as asterisk, convert to semicolon
    #     int = gsub(" \\* ", "\\:", basis.set[[i]][1])
    #     
    #     # Get all combinations of interactions
    #     all.ints = sapply(strsplit(int, ":"), function(x) { 
    #       
    #       datf = expand.grid(rep(list(x), length(x)), stringsAsFactors = FALSE)
    #       
    #       datf = datf[apply(datf, 1, function(x) !any(duplicated(x))), ]
    #       
    #       apply(datf, 1, function(x) paste(x, collapse = ":"))
    #       
    #     } )
    #     
    #     row.num = which(attr(terms(basis.mod.new), "term.labels") %in% all.ints) + 1
    #       
    #     }
    #   
    #   } else {
    #     
    #     row.num = which(basis.set[[i]][1] == basis.mod.new$varNames)
    #      
    #   }
    
    # Get row number of coefficient table
    row.num = length(attr(terms(basis.mod.new), "term.labels")) + 1
    
    # Return new coefficient table
    ret = if(any(class(basis.mod.new) %in% c("lmerMod", "merModLmerTest"))) {
      
      coef.table = suppressMessages(summary(basis.mod.new)$coefficients)
      
      # Get Kenward-Rogers approximation of denominator degrees of freedom
      kr.ddf = suppressMessages(get_ddf_Lb(basis.mod.new, fixef(basis.mod.new)))
      
      # Compute p-values based on t-distribution and ddf
      kr.p = 2 * (1 - pt(abs(coef.table[, "t value"]), kr.ddf))[row.num]
      
      # Combine with coefficients from regular ouput
      data.frame(
        t(coef.table[row.num, 1:2]),
        kr.ddf,
        coef.table[row.num, 3],
        kr.p,
        row.names = NULL
      )
      
    } else if(any(class(basis.mod.new) %in% c("lm", "glm", "negbin", "pgls", "glmerMod", "glmmadmb")))
      
      as.data.frame(t(unname(summary(basis.mod.new)$coefficients[row.num, ]))) else
      
        as.data.frame(t(unname(summary(basis.mod.new)$tTable[row.num, ])))  
    
    # Add df if summary table does not return
    if(length(ret) != 5 & any(class(basis.mod.new) %in% c("lm", "glm", "negbin", "pgls"))) 
      
      ret = cbind(ret[1:2], summary(basis.mod.new)$df[2], ret[3:4]) else
        
        if(length(ret) != 5 & any(class(basis.mod.new) %in% c("glmmadmb"))) 
          
          ret = cbind(ret[1:2], summary(basis.mod.new)$n, ret[3:4]) else
            
            if(length(ret) != 5)
              
              ret = cbind(ret[1:2], NA, ret[3:4])
      
    # Rename columns 
    names(ret) = c("estimate", "std.error", "df", "crit.value", "p.value")
    
    # Adjust p-value based on Shipley 2013
    if(adjust.p == TRUE) {
      
      if(any(class(basis.mod.new) %in% c("lme", "glmmPQL"))) {
        
        t.value = summary(basis.mod.new)$tTable[row.num, 4] 
        
        ret[5] = 2*(1 - pt(abs(t.value), nobs(basis.mod.new) - sum(apply(basis.mod.new$groups, 2, function(x) length(unique(x))))))
        
      } else if(any(class(basis.mod.new) %in% c("lmerMod", "glmerMod", "merModLmerTest"))) {
        
        z.value = coef.table[row.num, "t value"]
        
        ret[5] = 2*(1 - pt(abs(z.value), nobs(basis.mod.new) - sum(summary(basis.mod.new)$ngrps))) 
        
        } else if(any(class(basis.mod.new) %in% c("glmmadmb"))) {
          
          z.value = summary(basis.mod.new)$coefficients[row.num, 3]
          
          ret[5] = 2*(1 - pt(abs(z.value), nobs(basis.mod.new) - sum(summary(basis.mod.new)$npar))) 
      
        }
      
    }
  
    if(.progressBar == TRUE) setTxtProgressBar(pb, i)
    
    # Modify rhs if number of characters exceeds 20
    rhs.new = 
      
        if(length(basis.set[[i]]) < 3) rhs else {
          
          if(conditional == FALSE) 
            
            paste0(basis.set[[i]][1], " + ...") else 
              
              paste(basis.set[[i]][c(1, 3:length(basis.set[[i]]))], collapse = " + ")
          
        }
    
    # Bind in d-sep metadata
    data.frame(missing.path = paste(basis.set[[i]][2], " ~ ", rhs.new, sep = ""), ret)
    
  } ) ) else
    
    pvalues.df = data.frame(missing.path = NA, estimate = NA, std.error = NA, DF = NA, crit.value = NA, p.value = NA)
  
  # Identify duplicate tests from intermediate endogenous variables
  dup = names(basis.set)
  
  # Return lowest P-value
  pvalues.df = do.call(rbind, lapply(unique(dup), function(x) {
    
    if(length(dup[dup == x]) > 1) 
      
      warning("Some d-sep tests are non-symmetrical. The most conservative P-value has been returned. Stay tuned for future developments...")
    
    pvalues.df[as.numeric(x), ][which.min(pvalues.df[as.numeric(x), "p.value"]), ]
    
  }
  
  ) )
  
  # Set degrees of freedom as numeric
  pvalues.df$df = round(as.numeric(pvalues.df$df), 1)
  
  if(!is.null(pb)) close(pb)  
  
  if(any(grepl("...", pvalues.df$missing.path))) 
    
    message("Conditional variables have been omitted from output table for clarity (or use argument conditional = T)")
  
  rm(dup)
  
  return(pvalues.df)
  
}