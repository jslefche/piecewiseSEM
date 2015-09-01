get.random.formula = function(model, rhs, modelList, drop.terms = NULL) {
  
  if(class(rhs) == "formula") rhs = deparse(rhs)
  
  # Get random formula from model
  random.formula = if(any(class(model) %in% c("lme", "glmmPQL")))
    
    deparse(model$call$random) else
      
      if(any(class(model) %in% c("lmerMod", "merModLmerTest", "glmerMod")))
        
        paste("(", findbars(formula(model)), ")", collapse = " + ")
  
  # Get random structure(s)
  random.structure = if(any(class(model) %in% c("lme", "glmmPQL")))
    
    gsub(".*\\|(.*)?", "\\1", random.formula) else 
      
      if(any(class(model) %in% c("lmerMod", "merModLmerTest", "glmerMod")))
        
        sapply(strsplit(random.formula, ".\\+.")[[1]], function(x)
          
          gsub(".*\\|(.*)\\)", "\\1", x) 
          
          )
  
  # Get random slopes in the model list, otherwise return vector of terms to drop
  random.slopes = if(any(class(model) %in% c("lme", "glmmPQL", "glmerMod", "merModLmerTest"))) 
    
    if(is.null(drop.terms)) {
      
      sapply(1:length(modelList), function(j) {
        
        if(any(class(modelList[[j]]) %in% c("glmmPQL")))
          
          sapply(model$coefficients$random, function(k) colnames(k)[!colnames(k) %in% c("(Intercept)")])
        
        else if(any(class(modelList[[j]]) %in% c("lme", "glmerMod", "merModLmerTest")))
          
          colnames(ranef(modelList[[j]]))[!colnames(ranef(modelList[[j]])) %in% "(Intercept)"]
        
      } )
      
    } else drop.terms 

  random.slopes = unname(random.slopes[!duplicated(random.slopes)])
  
  # Define new random slopes 
  new.random.slopes = random.slopes[which(random.slopes %in% unlist(strsplit(rhs, ".\\+.")))]
  
  if(length(new.random.slopes) == 0) new.random.slopes = 1
  
  # Replace random slopes if any variables in model formula appear in random slopes  
  if(!length(random.slopes) == 0 & is.null(drop.terms)) {

    if(any(class(model) %in% c("lme", "glmmPQL")))
      
      paste("~ ", 
            paste(new.random.slopes, collapse = " + "),
            " | ",
            random.structure) 
    
    else if(any(class(model) %in% c("glmerMod", "merModLmerTest")))
      
      paste(
        sapply(random.structure, function(x)
          
          paste("(", paste(new.random.slopes, collapse = " + "), " | ", x, ")") ),
        
        collapse = " + ")
    
  } else if(!length(random.slopes) == 0 & !is.null(drop.terms)) {
    
    if(any(class(model) %in% c("lme", "glmmPQL")))
      
      paste("~ ",
            paste(new.random.slopes, collapse = " + "),
            " | ",
            random.structure) 
    
    else if(any(class(model) %in% c("glmerMod", "merModLmerTest")))
      
      paste(
      sapply(random.structure, function(x)
        
        paste("(", paste(new.random.slopes, collapse = " + "), " | ", x, ")") ),
      
      collapse = " + ")
    
  } else if(length(random.slopes) == 0)
    
    random.formula

}