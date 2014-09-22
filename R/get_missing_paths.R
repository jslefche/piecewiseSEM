get.missing.paths = function(modelList, data, corr.errors = NULL, add.vars = NULL, 
                             grouping.vars = NULL, top.level.vars = NULL,
                             adjust.p = FALSE, basis.set = NULL, disp.conditional = FALSE,
                             .progressBar = TRUE) {
  
  if(is.null(basis.set)) { 
    
    basis.set = get.basis.set(modelList, corr.errors, add.vars)
    
    basis.set = filter.exogenous(modelList, basis.set, corr.errors, add.vars) 
    
  }
  
  if(.progressBar == T) pb = txtProgressBar(min = 0, max = length(basis.set), style = 3) else pb = NULL
  
  pvalues.df = do.call(rbind, lapply(seq_along(basis.set), function(i) {
    
    basis.mod = modelList[[match(basis.set[[i]][2], unlist(lapply(modelList, function(j) as.character(formula(j)[2]))))]]
    
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
    
    if(is.null(random.formula)) basis.mod = update(basis.mod, fixed.formula, data = data) else
      if(all(class(basis.mod) %in% c("lme", "glmmPQL"))) 
        basis.mod = update(basis.mod, fixed = formula(fixed.formula), random = formula(random.formula), data = data) else
          basis.mod = update(basis.mod, formula = formula(paste(fixed.formula, "+", random.formula, sep="", collapse="+")), data = data)
    
    if(any(class(basis.mod) %in% "lmerMod")) basis.mod = as(basis.mod, "merModLmerTest") 
    
    ###
    
    if(all(class(basis.mod) %in% c("lmerMod", "merModLmerTest"))) {
      x = try(suppressMessages(suppressWarnings(summary(basis.mod))$coefficients[1,5]), silent = T)
      if(class(x) == "try-error") stop("lmerTest did not converge, no p-values to report. Consider specifying lmerControl") }
    
    ###
    
    if(!grepl(":", basis.set[[i]][1])) rowname = basis.set[[i]][1] else {
      int = attr(terms(basis.mod), "term.labels")[grepl(":", attr(terms(basis.mod), "term.labels"))]
      rowname = int[sapply(lapply(int, function(j) unlist(strsplit(j, ":"))), 
                           function(k) all(unlist(strsplit(basis.set[[i]][1], ":")) %in% k))]
    }
    
    if(adjust.p == TRUE) {
      if(all(class(basis.mod) %in% c("lm", "glm", "negbin"))) {
        p = summary(basis.mod)$coefficients[pmatch(rowname, rownames(summary(basis.mod)$coefficients)), 3] 
        df = basis.mod$df.residual
      } else if(all(class(basis.mod) %in% c("lme", "glmmPQL"))) {
        t.value = summary(basis.mod)$tTable[pmatch(rowname, rownames(summary(basis.mod)$tTable)), 4] 
        p = 2*(1 - pt(abs(t.value), nobs(basis.mod) - sum(apply(basis.mod$groups,2,function(x) length(unique(x))))))
        df = summary(basis.mod)$tTable[pmatch(rowname, rownames(summary(basis.mod)$tTable)), 3]
      } else if(all(class(basis.mod) %in% c("glmerMod"))) {
        z.value = summary(basis.mod)$coefficients[pmatch(rowname, rownames(summary(basis.mod)$coefficients)), 4]
        p = 2*(1 - pt(abs(z.value), nobs(basis.mod) - sum(summary(basis.mod)$ngrps)))
        df = "-"
      } else if(all(class(basis.mod) %in% c("merModLmerTest"))) {
        t.value = summary(basis.mod)$coefficients[pmatch(rowname, rownames(summary(basis.mod)$coefficients)), 4]
        p = 2*(1 - pt(abs(t.value), nobs(basis.mod) - sum(summary(basis.mod)$ngrps)))
        df = summary(basis.mod)$coefficients[pmatch(rowname, rownames(summary(basis.mod)$coefficients)), 3] }  
    } else if(adjust.p == FALSE) {
      if(all(class(basis.mod) %in% c("lm", "glm", "negbin"))) {
        p = summary(basis.mod)$coefficients[pmatch(rowname, rownames(summary(basis.mod)$coefficients)), 4]
        df = basis.mod$df.residual
      } else if(all(class(basis.mod) %in% c("lme", "glmmPQL"))) {
        p = summary(basis.mod)$tTable[pmatch(rowname, rownames(summary(basis.mod)$tTable)), 5]
        df = summary(basis.mod)$tTable[pmatch(rowname, rownames(summary(basis.mod)$tTable)), 3]
      } else if(class(basis.mod) %in% c("glmerMod")) {
        p = summary(basis.mod)$coefficients[pmatch(rowname, rownames(summary(basis.mod)$coefficients)), 4]
        df = "-"
      } else if(class(basis.mod) %in% c("merModLmerTest")) {
        p = summary(basis.mod)$coefficients[pmatch(rowname, rownames(summary(basis.mod)$coefficients)), 5]
        df = summary(basis.mod)$coefficients[pmatch(rowname, rownames(summary(basis.mod)$coefficients)), 3]  }
    }
    
    if(.progressBar == TRUE) setTxtProgressBar(pb, i)
    
    if(disp.conditional == FALSE)
      data.frame(
        missing.path = paste(basis.set[[i]][2], "<-", paste(basis.set[[i]][1], collapse = "+")), 
        df = df,
        p.value = round(p, 3))  else
          data.frame(
            missing.path = paste(basis.set[[i]][2], "<-", paste(basis.set[[i]][1], collapse = "+")), 
            df = df, 
            conditional.on = paste(basis.set[[i]][3:length(basis.set[[i]])], collapse = ","),
            p.value = round(p, 3))
    
  } ) )
  
  if(!is.null(pb)) close(pb)  
  
  return(pvalues.df)
  
}