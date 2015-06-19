predict.sem = function(modelList, newdata, se.fit = FALSE, ...) {
  
  # Send newdata to each model in the model list and return output as a data.frame 
  predict.df = do.call(cbind, lapply(modelList, function(i) {
    
    # Get model predictions
    if(any(class(i) %in% c("lm", "glm", "neg.bin", "gls", "pgls")) )
       
       predict.df = predict(i, newdata, se.fit = se.fit, ...) else
         
         if(any(class(i) %in% c("lme", "glmmPQL"))) 
           
           predict.df = predict(i, newdata, level = 0, ...) else
             
             if(any(class(i) %in% c("lmerMod", "glmerMod", "merModTest")))
               
               predict.df = predict(i, newdata, re.form = NA, ...)
    
             
  # If se.fit = TRUE for mixed models, calculate standard errors based on fixed-effects only
  if(se.fit == TRUE & any(class(i) %in% c("lme", "glmmPQL", "lmerMod", "glmerMod", "merModTest"))) {
    
    # Bind in predictions to new data
    newdata = data.frame(newdata, predict.df)
    
    colnames(newdata)[ncol(newdata)] = all.vars(formula(i))[1]
    
    if(any(class(i) %in% c("lme", "glmmPQL"))) {
      
      Dmat = model.matrix(formula(i)[-2], newdata) 
      
      pvar = sqrt(diag(Dmat %*% vcov(i) %*% t(Dmat)))
      
    } else {
      
      Dmat = model.matrix(terms(i), newdata)
      
      pvar = diag(Dmat %*% tcrossprod(vcov(i), Dmat))
      
    }
    
    # Return list with predicted errors
    predict.df = list(fit = predict.df, se.fit = pvar)
    
  }
      
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