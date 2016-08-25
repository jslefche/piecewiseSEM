sem.missing.paths = function(
 
  modelList, data, conditional = FALSE, corr.errors = NULL, add.vars = NULL, grouping.vars = NULL, 
  adjust.p = FALSE, basis.set = NULL, model.control = NULL, .progressBar = TRUE
  
  ) {
  
  # Get basis set
  if(is.null(basis.set)) basis.set = suppressWarnings(sem.basis.set(modelList, corr.errors, add.vars))

  # Add progress bar
  if(.progressBar == T & length(basis.set) > 0) pb = txtProgressBar(min = 0, max = length(basis.set), style = 3) else pb = NULL
  
  # Perform d-sep tests
  if(length(basis.set) > 0) pvalues.df = do.call(rbind, lapply(1:length(basis.set), function(i) {
    
    # Get basis model from which to build the d-sep test
    basis.mod = modelList[[which(sapply(modelList, function(j) {
      
      if(any(class(j) == "pgls")) j = j$formula
      
      gsub(" ", "", rownames(attr(terms(j), "factors"))[1]) == basis.set[[i]][2]
    
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
    basis.mod.new = suppressMessages(suppressWarnings(
      
      if(is.null(random.formula) | class(basis.mod) == "glmmadmb") 
      
        update(basis.mod, formula(paste(basis.set[[i]][2], " ~ ", rhs)), control = control, data = data) else
        
          if(any(class(basis.mod) %in% c("lme", "glmmPQL"))) 
          
            update(basis.mod, fixed = formula(paste(basis.set[[i]][2], " ~ ", rhs)), random = random.formula, control = control, data = data) else
            
              update(basis.mod, formula = formula(paste(basis.set[[i]][2], " ~ ", rhs, " + ", random.formula)), control = control, data = data) 
      
    ) )

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
    
    # Return new coefficient table
    ret = if(any(class(basis.mod.new) %in% c("lmerMod", "merModLmerTest"))) {
      
      coef.table = suppressMessages(summary(basis.mod.new)$coefficients)

      # Get P-values baesd on Kenward-Rogers approximation of denominator degrees of freedom
      basis.mod.drop = update(basis.mod.new, as.formula(paste("~ . -", basis.set[[i]][1])))
      
      kr.p = KRmodcomp(basis.mod.new, basis.mod.drop)
      
      # Combine with coefficients from regular ouput
      data.frame(
        t(coef.table[nrow(coef.table), 1:2]),
        kr.p$test$ddf[1],
        coef.table[nrow(coef.table), 3],
        kr.p$test$p.value[1],
        row.names = NULL
      )
      
    } else if(any(class(basis.mod.new) %in% c("lm", "glm", "negbin", "pgls", "glmerMod", "glmmadmb"))) {
      
      coef.table = summary(basis.mod.new)$coefficients
      
      as.data.frame(t(unname(coef.table[nrow(coef.table), ]))) 
      
      } else {
      
        coef.table = summary(basis.mod.new)$tTable
        
        as.data.frame(t(unname(coef.table[nrow(coef.table), ])))  
        
      }
    
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
        
        t.value = coef.table[nrow(coef.table), 4] 
        
        ret[5] = 2*(1 - pt(abs(t.value), nobs(basis.mod.new) - sum(apply(basis.mod.new$groups, 2, function(x) length(unique(x))))))
        
      } else if(any(class(basis.mod.new) %in% c("lmerMod", "glmerMod"))) {
        
        z.value = coef.table[nrow(coef.table), 3]
        
        ret[5] = 2*(1 - pt(abs(z.value), nobs(basis.mod.new) - sum(summary(basis.mod.new)$ngrps))) 
        
        } else if(any(class(basis.mod.new) %in% c("glmmadmb"))) {
          
          z.value = coef.table[nrow(coef.table), 3]
          
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
  
  if(any(grepl("...", pvalues.df$missing.path)) & conditional != TRUE) 
    
    message("Conditional variables have been omitted from output table for clarity (or use argument conditional = T)")
  
  # rm(dup)
  
  # Assign significance indicators
  sig = sapply(pvalues.df$p.value, function(y) {
    
    ifelse(y > 0.01 & y < 0.05, "*", 
           ifelse(y > 0.001 & y <= 0.01, "**",
                  ifelse(y <= 0.001, "***", "")
           )
    )
    
  } )
  
  pvalues.df = cbind(pvalues.df, sig)
  
  colnames(pvalues.df)[ncol(pvalues.df)] = ""
  
  return(pvalues.df)
  
}