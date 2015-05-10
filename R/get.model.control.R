get.model.control = function(model, model.control) {
  
  if(length(model.control)>5) model.control = list(model.control)
  
  # Match model control list with appropriate model class for basis model
  if(is.null(model.control)) {
    
    if(any(class(model) %in% "glm")) glm.control() else
      
      if(any(class(model) %in% "gls")) gls.control() else
      
        if(class(model) %in% c("lme", "glmmPQL")) lmeControl() else 
        
          if(class(model) %in% c("lmerMod", "merModLmerTest")) lmerControl() else
            
            if(class(model) %in% c("glmerMod")) glmerControl()
          
  } else {
    
    classes = lapply(model.control, function(i) gsub("(.*)Control", "\\1", class(i))[1] )
    
    if(any(classes %in% gsub("(.*)Mod", "\\1", class(model))))
      
      model.control[[which(classes %in% gsub("(.*)Mod", "\\1", class(model)))]] else
        
        if(class(model) == "glm")
          
          model.control[[which(sapply(model.control, length) == 3)]] else
            
            if(class(model) == "gls")
              
              model.control[[which(sapply(model.control, length) == 13)]] else
                
                if(class(model) %in% c("lme", "glmmPQL"))
                  
                  model.control[[which(sapply(model.control, length) == 15)]]
            
  }

}