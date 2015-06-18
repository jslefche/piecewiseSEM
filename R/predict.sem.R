predict.sem = function(modelList, newdata, ...) {
  
  # Send newdata to each model in the model list and return output as a data.frame 
  predict.df = do.call(cbind, lapply(modelList, function(i) {
    
    # Get model predictions
    if(any(class(i) %in% c("lm", "glm", "neg.bin", "gls", "pgls")) )
       
       predict.df = predict(i, newdata, ...) else
         
         if(any(class(i) %in% c("lme", "glmmPQL"))) 
           
           predict.df = predict(i, newdata, level = 0, ...) else
             
             if(any(class(i) %in% c("lmerMod", "glmerMod", "merModTest")))
               
               predict.df = predict(i, newdata, re.form = NA, ...)
    
    # If predictions are stored in a list, bind columns
    if(class(predict.df) == "list") predict.df = do.call(data.frame, predict.df[1:2]) else
      
      predict.df = data.frame(predict.df)
  
    # Name columns
    if(ncol(predict.df) == 1)
      
      colnames(predict.df) = paste(all.vars(formula(i))[1], "fit", sep = ".") else 
        
        colnames(predict.df) = paste(all.vars(formula(i))[1], colnames(predict.df), sep = ".")
    
    # Return predictions
    return(predict.df)
      
  } )
  
  )
  
  # Bind in newdata
  cbind(newdata, predict.df)
  
}