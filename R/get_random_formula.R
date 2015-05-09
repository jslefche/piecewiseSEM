get.random.formula = function(model, rhs, modelList, drop.terms = NULL) {
  
  if(class(rhs) == "formula") rhs = deparse(rhs)
  
  # Get random formula from model
  random.formula = if(any(class(model) %in% c("lme", "glmmPQL")))
    
    model$call$random else
      
      if(any(class(model) %in% c("lmerMod", "merModLmerTest", "glmerMod")))
        
        paste("(", findbars(formula(model)), ")", collapse = " + ")
  
  # Get random structure(s)
  random.structure = if(any(class(model) %in% c("lme", "glmmPQL")))
    
    gsub(".*\\|(.*)?", "\\1", random.formula)[2] else 
      
      if(any(class(model) %in% c("lmerMod", "merModLmerTest", "glmerMod")))
        
        sapply(strsplit(deparse(random.formula), ".\\+.")[[1]], function(x)
          
          gsub(".*\\|(.*)\\)", "\\1", x) 
          
          )
  
  # Get random slopes in the model list, otherwise return vector of terms to drop
  random.slopes = if(any(class(model) %in% c("lme", "glmmPQL", "glmerMod", "merModLmerTest"))) 
    
    if(is.null(drop.terms)) 
      
      suppressWarnings(unlist(lapply(
        
        modelList[sapply(modelList, class) %in% c("lme", "glmmPQL", "glmerMod", "merModLmerTest")], 
        
        function(i) 
          
          colnames(do.call(cbind, ranef(i)))[!colnames(do.call(cbind, ranef(i))) %in% "(Intercept)"]
        
        ) ) ) else drop.terms 

  random.slopes = random.slopes[!duplicated(random.slopes)]
  
  # Define new random slopes 
  new.random.slopes = random.slopes[which(random.slopes %in% unlist(strsplit(rhs, ".\\+.")))]
  if(length(new.random.slopes) == 0) new.random.slopes = 1
  
  # Replace random slopes if any variables in model formula appear in random slopes  
  if(!length(random.slopes) == 0 & is.null(drop.terms)) {

    if(class(model) %in% c("lme", "glmmPQL"))
      
      paste("~ ", 
            paste(new.random.slopes, collapse = " + "),
            " | ",
            random.structure) 
    
    else if(class(model) %in% c("glmerMod", "merModLmerTest"))
      
      paste(
        sapply(random.structure, function(x)
          
          paste("(", paste(new.random.slopes, collapse = " + "), " | ", x, ")") ),
        
        collapse = " + ")
    
  } else if(!length(random.slopes) == 0 & !is.null(drop.terms)) {
    
    if(class(model) %in% c("lme", "glmmPQL"))
      
      paste("~ ",
            paste(new.random.slopes, collapse = " + "),
            " | ",
            random.structure) 
    
    else if(class(model) %in% c("glmerMod", "merModLmerTest"))
      
      paste(
      sapply(random.structure, function(x)
        
        paste("(", paste(new.random.slopes, collapse = " + "), " | ", x, ")") ),
      
      collapse = " + ")
    
  } else if(length(random.slopes) == 0)
    
    random.formula

}