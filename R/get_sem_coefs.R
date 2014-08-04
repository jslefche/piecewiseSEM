get.sem.coefs = function(modelList, standardized = FALSE) {

  names(modelList) = NULL
 
  if(standardized == FALSE) {
    df = do.call(rbind, lapply(modelList, function(i) {
      
      if(all(class(i) %in% c("lm", "glm", "negbin", "glmerMod"))) {
        tab = summary(i)$coefficients
        data.frame(path = paste(Reduce(paste, deparse(formula(i)[[2]])), "<-", rownames(tab)[-1]),
                   estimate = round(tab[-1, 1], 3),
                   std.error = round(tab[-1, 2], 3),
                   p.value = round(tab[-1, 4], 3), 
                   row.names = NULL)
        
      } else if(all(class(i) %in% c("lme", "glmmPQL"))) {
        tab = summary(i)$tTable
        data.frame(path = paste(Reduce(paste, deparse(formula(i)[[2]])), "<-", rownames(tab)[-1]),
                   estimate = round(tab[-1, 1], 3),
                   std.error = round(tab[-1, 2], 3),
                   p.value = round(tab[-1, 5], 3), 
                   row.names = NULL)
        
      } else if(all(class(i) %in% c("lmerMod", "merModLmerTest"))) {
        tab = summary(as(i, "merModLmerTest"))$coefficients
        data.frame(path = paste(Reduce(paste, deparse(formula(i)[[2]])), "<-", rownames(tab)[-1]),
                   estimate = round(tab[-1, 1], 3),
                   std.error = round(tab[-1, 2], 3),
                   p.value = round(tab[-1, 5], 3), 
                   row.names = NULL) } } ) )
    
  } else if(standardized == TRUE) {
    
    df = do.call(rbind, lapply(modelList, function(i) {
      
      if(all(class(i) %in% c("lm", "glm", "negbin"))) model.data = i$model else 
        if(all(class(i) %in% c("lme", "glmmPQL"))) model.data = i$data else
          if(all(class(i) %in% c("lmerMod", "merModLmerTest", "glmerMod"))) model.data = i@frame

      vars.toscale = if(all(class(i) %in% c("lm", "lme"))) rownames(attr(i$terms, "factors")) else
        if(any(class(i) %in% c("glm", "negbin", "glmmPQL"))) {
          message("Model is not gaussian: keeping response on original scale")
          rownames(attr(i$terms, "factors"))[-1] } else 
            if(all(class(i) %in% c("merModLmerTest")))
              rownames(attr(attr(i@frame, "terms"), "factors"))[!rownames(attr(attr(i@frame, "terms"), "factors")) %in% names(i@flist)] else
                if(all(class(i) %in% c("glmerMod"))) {
                  message("Model is not gaussian: keeping response on original scale")
                  rownames(attr(attr(i@frame, "terms"), "factors"))[!rownames(attr(attr(i@frame, "terms"), "factors")) %in% names(i@flist)][-1] }

      model.data[, vars.toscale] = apply(model.data[, vars.toscale, drop = F], 2, function(x) if(is.numeric(x)) scale(x) else x )
      
      model = update(i, data = model.data)
      
      if(class(model) == "lmerMod") model = as(model, "merModLmerTest")
      
      if(all(class(model) %in% c("lm", "glm", "negbin", "glmerMod"))) {
        tab = summary(model)$coefficients
        data.frame(path = paste(Reduce(paste, deparse(formula(model)[[2]])), "<-", rownames(tab)[-1]),
                   estimate = round(tab[-1, 1], 3),
                   std.error = round(tab[-1, 2], 3),
                   p.value = round(tab[-1, 4], 3), 
                   row.names = NULL)
        
      } else if(all(class(model) %in% c("lme", "glmmPQL"))) {
        tab = summary(model)$tTable
        data.frame(path = paste(Reduce(paste, deparse(formula(model)[[2]])), "<-", rownames(tab)[-1]),
                   estimate = round(tab[-1, 1], 3),
                   std.error = round(tab[-1, 2], 3),
                   p.value = round(tab[-1, 5], 3), 
                   row.names = NULL)
        
      } else if(all(class(model) %in% c("merModLmerTest"))) {
        tab = summary(model)$coefficients
        data.frame(path = paste(Reduce(paste, deparse(formula(model)[[2]])), "<-", rownames(tab)[-1]),
                   estimate = round(tab[-1, 1], 3),
                   std.error = round(tab[-1, 2], 3),
                   p.value = round(tab[-1, 5], 3), 
                   row.names = NULL) } 
    } ) )
  }
  
  df[order(df$p.value),]
}