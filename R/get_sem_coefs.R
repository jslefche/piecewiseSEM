get.sem.coefs = function(modelList, standardized = FALSE, corr.errors = NULL) {
  
  names(modelList) = NULL
  
  if(standardized == FALSE) {
    dataframe = do.call(rbind, lapply(modelList, function(i) {
      
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
    
    dataframe = do.call(rbind, lapply(modelList, function(i) {
      
      vars.to.scale = if(all(class(i) %in% c("lm", "lme"))) rownames(attr(i$terms, "factors")) else
        if(any(class(i) %in% c("glm", "negbin", "glmmPQL"))) {
          message("Model is not gaussian: keeping response on original scale")
          rownames(attr(i$terms, "factors"))[-1] } else 
            if(all(class(i) %in% c("merModLmerTest")))
              rownames(attr(attr(i@frame, "terms"), "factors"))[!rownames(attr(attr(i@frame, "terms"), "factors")) %in% names(i@flist)] else
                if(all(class(i) %in% c("glmerMod"))) {
                  message("Reponse is not modeled to a gaussian distribution: keeping response on original scale")
                  rownames(attr(attr(i@frame, "terms"), "factors"))[!rownames(attr(attr(i@frame, "terms"), "factors")) %in% names(i@flist)][-1] }
      
      form = gsub(" ", "", unlist(strsplit(paste(format(formula(i)), collapse = ""), "\\+|\\~")))
 
      if(any(grepl("\\*|\\:", form))) {
        ints = strsplit(form[grepl("\\*|\\:", form)], "\\*")
        ints.scaled = lapply(ints, function(j) paste(paste("scale(",j,")"), collapse = "*"))
        form[grepl("\\*|\\:", form)] = unlist(ints.scaled)
        vars.to.scale=vars.to.scale[-match(unique(unlist(ints)),vars.to.scale)] }
        
      if(length(vars.to.scale) > 0) form[match(vars.to.scale, form)] = paste("scale(", form[match(vars.to.scale, form)], ")")
      
      form = if(length(form) > 2)
        paste(paste(form[1], "~", form[2], "+"), paste(form[-c(1:2)], collapse = "+")) else
          paste(paste(form[1], "~", form[2]))
      
      if(all(class(i) %in% c("lme", "glmmPQL")))
        model = update(i, fixed = formula(form)) else
          model = update(i, formula = form)
      
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
  
  dataframe = dataframe[order(dataframe$p.value),]
  
  if(!is.null(corr.errors)) {
    dataframe = rbind(dataframe, do.call(rbind, lapply(corr.errors, function(j) {
      
      corr.vars = gsub(" ", "", unlist(strsplit(j,"~~")))
      
      if(all(corr.vars %in% unlist(lapply(modelList, function(i) as.character(formula(i)[2]))))) {
        
        resid.data = get.partial.resid(y = corr.vars[1], x = corr.vars[2], modelList)
        
        data.frame(
          path = j,
          estimate = round(cor(resid.data)[1,2], 3),
          std.error = NA,
          p.value = round(1 - pt((cor(resid.data)[1, 2] * sqrt(nrow(resid.data) - 2))/(sqrt(1 - cor(resid.data)[1, 2]^2)),nrow(resid.data)-2), 3),
          row.names = NULL) 
        
      } else {
        
        model.data = do.call(cbind, lapply(modelList, function(i)           
          if(all(class(i) %in% c("lm", "glm", "negbin"))) model.data = i$model else 
            if(all(class(i) %in% c("lme", "glmmPQL"))) model.data = i$data else
              if(all(class(i) %in% c("lmerMod", "merModLmerTest", "glmerMod"))) model.data = i@frame ) )
        
        model.data = model.data[!duplicated(colnames(model.data))]
        
        data.frame(
          path = j,
          estimate = round(cor(model.data[, corr.vars[1]], model.data[, corr.vars[2]]), 3),
          std.error = "-",
          p.value = round(1 - pt((cor(model.data[, corr.vars[1]], model.data[, corr.vars[2]]) * sqrt(nrow(model.data) - 2))/(sqrt(1 - cor(model.data[, corr.vars[1]], model.data[, corr.vars[2]])^2)),nrow(model.data)-2), 3),
          row.names = NULL)
        
      }
      
    } ) ) ) }
  
  return(dataframe)
  
}