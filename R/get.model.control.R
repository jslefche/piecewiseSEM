# Fixes issues with nlme updated; thanks to Martin Maechler, 2016-01-19

get.model.control = function(model, model.control) {
  
  model.class = if(inherits(model, "merModLmerTest")) "lmerMod" else class(model)
  
  # Match model control list with appropriate model class for basis model
  if(is.null(model.control)) {
    
    if(inherits(model, "glm")) glm.control()
    
    else if(inherits(model, "gls")) glsControl()
    
    else if(any(class(model) %in% c("lme", "glmmPQL"))) lmeControl()
    
    else if(any(class(model) %in% c("lmerMod", "merModLmerTest"))) lmerControl()
    
    else if(inherits(model, "glmerMod")) glmerControl()
    
  } else {
    
    ## assume model.control = a *list* of control lists
    control.classes = lapply(model.control, function(i)
      
      gsub("(.*)Control", "\\1", class(i))[1] )
    
    if(any(control.classes %in% (M.cl <- gsub("(.*)Mod", "\\1", model.class))))
      
      model.control[[control.classes %in% M.cl]] else
        
        if(inherits(model, "glm"))
          
          model.control[[sapply(model.control, length) >= 3]] else
            
            if(inherits(model, "gls")) 
              
              model.control[[sapply(model.control, length) >= 13]] else
                
                if(any(class(model) %in% c("lme", "glmmPQL"))) 
                  
                  model.control[[sapply(model.control, length) >= 15]]
    
  }
  
}