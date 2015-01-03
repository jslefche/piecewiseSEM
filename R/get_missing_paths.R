get.missing.paths = function(modelList, data, corr.errors = NULL, add.vars = NULL, 
                             grouping.vars = NULL, top.level.vars = NULL,
                             adjust.p = FALSE, basis.set = NULL, disp.conditional = FALSE,
                             model.control = NULL, .progressBar = TRUE) {
  
  if(is.null(basis.set)) { 
    
    basis.set = get.basis.set(modelList, corr.errors, add.vars)
    
    basis.set = filter.exogenous(modelList, basis.set, corr.errors, add.vars) 
    
  }
  
  if(.progressBar == T) pb = txtProgressBar(min = 0, max = length(basis.set), style = 3) else pb = NULL
  
  pvalues.df = do.call(rbind, lapply(seq_along(basis.set), function(i) {
    
    basis.mod = modelList[[match(basis.set[[i]][2], unlist(lapply(modelList, function(j) as.character(formula(j)[2]))))]]
    
    #### Need to fix getting random effects structure    
    
    fixed.formula = paste(basis.set[[i]][2], "~", paste(basis.set[[i]][c(1, 3:length(basis.set[[i]]))], collapse = "+"))
    
    random.formula = if(all(class(basis.mod) %in% c("lme", "glmmPQL"))) 
      gsub("^.*\\|(.*?)", "\\1", basis.mod$call[[4]])[2] else
        if(all(class(basis.mod) %in% c("lmerMod", "merModLmerTest", "glmerMod")))
          gsub("^.*\\|(.*?)", "\\1", unlist(lapply(findbars(formula(basis.mod)),format))) else 
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
    
    if(!is.null(grouping.vars) & 
         any(top.level.vars %in% gsub(".*\\((.*)\\)", "\\1", unlist(as.character(formula(basis.mod)[2])))))
      data = suppressWarnings(aggregate(data, by = lapply(grouping.vars, function(i) data[ ,i]), mean, na.rm = T))
    
    if(length(model.control)>10) model.control = list(model.control)
    
    if(is.null(model.control)) {
      if(class(basis.mod) %in% c("lme", "glmmPQL")) control = lmeControl() else 
        if(class(basis.mod) %in% c("lmerMod", "merModLmerTest")) control = lmerControl() else
          if(class(basis.mod) %in% c("glmerMod")) control = glmerControl()} else {
            if(!is.null(model.control) & class(basis.mod) %in% c("lme", "glmmPQL"))
              control = model.control[[which(sapply(lapply(model.control, function(x) attr(x, "class")), is.null))]] else
                if(!is.null(model.control) & class(basis.mod) %in% c("lmerMod", "merModLmerTest"))
                  control = model.control[[which(sapply(lapply(model.control, function(x) attr(x, "class")), function(x) any(x %in% "lmerControl")))]] else
                    if(!is.null(model.control) & class(basis.mod) %in% c("glmerMod"))
                      control = model.control[[which(sapply(lapply(model.control, function(x) attr(x, "class")), function(x) any(x %in% "glmerControl")))]] }
    
    if(is.null(random.formula)) basis.mod = update(basis.mod, fixed.formula, data = data) else
      if(all(class(basis.mod) %in% c("lme", "glmmPQL"))) 
        basis.mod = update(basis.mod, fixed = formula(fixed.formula), random = formula(random.formula), control = control, data = data) else
          basis.mod = update(basis.mod, formula = formula(paste(fixed.formula, "+", random.formula, sep="", collapse="+")), control = control, data = data)
    
    if(any(class(basis.mod) %in% "lmerMod")) basis.mod = as(basis.mod, "merModLmerTest") 
    
    ###
    
    if(all(class(basis.mod) %in% c("lmerMod", "merModLmerTest"))) {
      x = try(suppressMessages(suppressWarnings(summary(basis.mod))$coefficients[1,5]), silent = T)
      if(class(x) == "try-error") stop("lmerTest did not converge, no p-values to report. Consider specifying model.control") }
    
    ###
    
    if(!grepl(":|\\*", basis.set[[i]][1])) row.num = which(basis.set[[i]][1] == attr(terms(basis.mod), "term.labels")) + 1 else
      row.num = which(grepl(":|\\*", attr(terms(basis.mod), "term.labels"))) + 1 
    
    ret = if(all(class(basis.mod) %in% c("lm", "glm", "negbin", "glmerMod", "merModLmerTest"))) 
      as.data.frame(t(unname(summary(basis.mod)$coefficients[row.num, ]))) else
        as.data.frame(t(unname(summary(basis.mod)$tTable[row.num, ])))
   
    if(length(ret) != 5) ret = cbind(ret[1:2], NA, ret[3:4])
    
    names(ret)=c("estimate","std.error","DF","crit.value","p.value")
    
    if(adjust.p == TRUE) {
      if(all(class(basis.mod) %in% c("lme", "glmmPQL"))) {
        t.value = summary(basis.mod)$tTable[row.num, 4] 
        ret[5] = 2*(1 - pt(abs(t.value), nobs(basis.mod) - sum(apply(basis.mod$groups,2,function(x) length(unique(x))))))
      } else if(all(class(basis.mod) %in% c("glmerMod", "merModLmerTest"))) {
        z.value = summary(basis.mod)$coefficients[row.num, 4]
        ret[5] = 2*(1 - pt(abs(z.value), nobs(basis.mod) - sum(summary(basis.mod)$ngrps))) } 
    } 
    
    if(disp.conditional == TRUE) 
      ret = cbind(missing.path = paste(basis.set[[i]][2], "<-", paste(basis.set[[i]][1], collapse = "+")), 
                  conditional.on = paste(basis.set[[i]][3:length(basis.set[[i]])], collapse = ","),
                  ret) else
                    ret = cbind(missing.path = paste(basis.set[[i]][2], "<-", paste(basis.set[[i]][1], collapse = "+")), 
                                ret)
    
    if(.progressBar == TRUE) setTxtProgressBar(pb, i)
    
    return(ret)
    
  } ) )
  
  if(!is.null(pb)) close(pb)  
  
  return(pvalues.df)
  
}