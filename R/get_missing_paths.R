get.missing.paths = function(modelList, adjust.p = FALSE, .progressBar = FALSE, basis.set = NULL, add.vars = NULL) {
  
  if(is.null(basis.set)) { 
    
    basis.set = get.basis.set(modelList, add.vars)
    
    basis.set = filter.exogenous(modelList, basis.set, add.vars) 
    
  }
  
  if(.progressBar == T & length(basis.set) > 1) pb = txtProgressBar(min = 0, max = length(basis.set), style = 3) else pb = NULL
  
  pvalues.df = do.call(rbind, lapply(seq_along(basis.set), function(i) {
    
    basis.mod = modelList[unlist(lapply(modelList, function(j) grepl(basis.set[[i]][2], formula(j)[2])))][[1]]
    
    fixed.formula = paste(basis.set[[i]][2], "~", paste(basis.set[[i]][c(1, 3:length(basis.set[[i]]))], collapse = "+"))
    
    random.formula = if(all(class(basis.mod) %in% c("lme", "glmmPQL"))) 
      gsub("^.*\\|(.*?)", "\\1", basis.mod$call[[4]])[2] else
        if(all(class(basis.mod) %in% c("lmerMod", "merModLmerTest", "glmerMod")))
          unlist(regmatches(format(formula(basis.mod)), gregexpr("(?<=\\|).*?(?=\\))", format(formula(basis.mod)), perl=TRUE))) else 
            NULL
    
    modelList.random.slopes = unlist(
      lapply(modelList, function(j) if(all(class(j) %in% c("lme", "glmmPQL"))) 
        gsub("\\|.*", "\\1", summary(j)$call[[4]][-1]) else
          if(all(class(j) %in% c("lmerMod", "merModLmerTest", "glmerMod"))) 
            gsub(".*\\((.*)\\|.*", "\\1", paste(formula(j)[[3]], collapse="")) else   
              NULL ) )
    
    random.slopes = unlist(lapply(basis.set[[i]][c(1, 3:length(basis.set[[i]]))], function(k) 
      if(any(grepl(k, modelList.random.slopes))) k else NULL) )
    if(is.null(random.slopes)) random.slopes = 1 else random.slopes = random.slopes
    
    random.formula = if(all(class(basis.mod) %in% c("lme", "glmmPQL")))
      paste(paste("~", paste(random.slopes, collapse="+"), "|", random.formula, sep = ""), collapse = "+") else
        if(all(class(basis.mod) %in% c("lmerMod", "merModLmerTest", "glmerMod")))
          paste(paste("(", paste(random.slopes, collapse="+"), "|", random.formula, ")", sep=""), collapse="+") else
            NULL
    
    if(is.null(random.formula)) basis.mod = update(basis.mod, fixed.formula) else
      if(all(class(basis.mod) %in% c("lme", "glmmPQL"))) 
        basis.mod = update(basis.mod, fixed = formula(fixed.formula), random = formula(random.formula)) else
          basis.mod = update(basis.mod, formula = formula(paste(fixed.formula, "+", random.formula, sep="", collapse="+")))
    
    if(class(basis.mod) %in% "lmerMod") basis.mod = as(basis.mod, "merModLmerTest") 
    
    if(!grepl(":", basis.set[[i]][1])) rowname = basis.set[[i]][1] else {
      int = attr(terms(basis.mod), "term.labels")[grepl(":", attr(terms(basis.mod), "term.labels"))]
      rowname = int[sapply(lapply(int, function(j) unlist(strsplit(j, ":"))), function(k) all(unlist(strsplit(basis.set[[i]][1], ":")) %in% k))]
    }
    
    if(adjust.p == TRUE) {
      if(all(class(basis.mod) %in% c("lm", "glm", "negbin"))) {
        p = summary(basis.mod)$coefficients[rowname,3] 
      } else if(all(class(basis.mod) %in% c("lme", "glmmPQL"))) {
        t.value = summary(basis.mod)$tTable[rowname,4] 
        p = 2*(1 - pt(abs(t.value), nobs(basis.mod) - sum(apply(basis.mod$groups,2,function(x) length(unique(x))))))
      } else if(all(class(basis.mod) %in% c("glmerMod"))) {
        z.value = summary(basis.mod)$coefficients[rowname,4]
        p = 2*(1 - pt(abs(z.value), nobs(basis.mod) - sum(summary(basis.mod)$ngrps)))
      } else if(all(class(basis.mod) %in% c("merModLmerTest"))) {
        t.value = summary(basis.mod)$coefficients[rowname,4]
        p = 2*(1 - pt(abs(t.value), nobs(basis.mod) - sum(summary(basis.mod)$ngrps))) }  
    } else if(adjust.p == FALSE) {
      if(all(class(basis.mod) %in% c("lm", "glm", "negbin", "glmerMod"))) {
        p = summary(basis.mod)$coefficients[rowname,4] 
      } else if(all(class(basis.mod) %in% c("lme", "glmmPQL"))) {
        p = summary(basis.mod)$tTable[rowname,5]
      } else if(class(basis.mod) %in% c("merModLmerTest")) {
        p = summary(basis.mod)$coefficients[rowname,5] }
    }
    
    if(.progressBar == TRUE) setTxtProgressBar(pb, i)
    
    data.frame(
      missing.path = paste(basis.set[[i]][2], "<-", paste(basis.set[[i]][1], collapse = "+")), 
      conditional.on = paste(basis.set[[i]][3:length(basis.set[[i]])], collapse = ","),
      p.value = round(p, 3))
    
  } ) )
  
  if(!is.null(pb)) close(pb)  
  
  return(pvalues.df)
  
}