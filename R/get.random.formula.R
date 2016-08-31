get.random.formula = function(model, rhs, modelList, dropterms = NULL) {
  
  if(class(rhs) == "formula") rhs = Reduce(paste, deparse(rhs))
  
  # Get random formula from model
  random.formula = if(any(class(model) %in% c("lme", "glmmPQL", "glmmadmb")))
    
    deparse(model$call$random) else
      
      if(any(class(model) %in% c("lmerMod", "merModLmerTest", "glmerMod", "glmmTMB")))
        
        paste("(", findbars(formula(model)), ")", collapse = " + ")
  
  # Get random structure(s)
  random.structure = if(any(class(model) %in% c("lme", "glmmPQL", "glmmadmb"))) {
    
    # If crossed random effects, extract random effects and store in a list
    if(grepl("list\\(", as.character(random.formula))) 
      
      as.list(
        gsub(
          "(.*)=.*",
          "\\1",
          strsplit(gsub("list\\((.*)\\)", "\\1", random.formula), ",")[[1]]
          )
        
      ) else gsub(".*\\|(.*)?", "\\1", random.formula) 
    
    } else 
      
      if(any(class(model) %in% c("lmerMod", "merModLmerTest", "glmerMod", "glmmTMB")))
        
        sapply(strsplit(random.formula, ".\\+.")[[1]], function(x)
          
          gsub(".*\\|(.*)\\)", "\\1", x) 
          
          )
  
  # Get random slopes in the model list, otherwise return vector of terms to drop
  random.slopes = if(any(class(model) %in% c("lme", "glmmPQL", "glmerMod", "merModLmerTest", "glmmadmb", "glmmTMB"))) 
    
    if(is.null(dropterms)) {
      
      unlist(lapply(1:length(modelList), function(j) {
        
        if(any(class(modelList[[j]]) %in% c("glmmPQL")))
          
          unlist(lapply(modelList[[j]]$coefficients$random, function(k) colnames(k)))
        
        else if(any(class(modelList[[j]]) %in% c("lme", "glmerMod", "merModLmerTest", "glmmadmb", "glmmTMB"))) 
          
          unlist(lapply(ranef(modelList[[j]]), function(k) colnames(k)))
        
      } ) )
      
    } else dropterms 

  random.slopes = unname(random.slopes[!duplicated(random.slopes) & random.slopes != "(Intercept)"])
  
  # Define new random slopes 
  new.random.slopes = random.slopes[which(random.slopes %in% unlist(strsplit(rhs, ".\\+.")))]
  
  if(length(new.random.slopes) == 0) new.random.slopes = 1 else new.random.slopes = paste0(new.random.slopes, collapse = " + ")
  
  # Replace random slopes if any variables in model formula appear in random slopes  
  if(length(random.slopes) != 0) {

    if(any(class(model) %in% c("lme", "glmmPQL", "glmmadmb")))
      
      if(is.list(random.structure)) {
        
        lapply(random.structure, function(x) formula(gsub("(.*)\\|", paste("~", new.random.slopes, "|"), x)) ) 
        
        } else {
        
        formula(
          paste("~ ", 
                new.random.slopes,
                " | ",
                random.structure) 
        )
    
    } else if(any(class(model) %in% c("glmerMod", "merModLmerTest", "glmmTMB")))
      
      # formula(
        paste(
          sapply(random.structure, function(x)
            paste("(", new.random.slopes, " | ", x, ")") ),
          collapse = " + ")
      # )
    
  } else if(length(random.slopes) == 0) {
      
      if(is.list(random.structure)) {
        
        lapply(random.structure, function(x) formula(gsub("(.*)\\|", paste("~", new.random.slopes, "|"), x)) ) 
        
      } else if(any(class(model) %in% c("lme", "glmmPQL"))) formula(random.formula)
    
    else random.formula
    
  }

}