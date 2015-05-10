partial.resid = function(.formula = y ~ x, modelList, model.control = NULL, plotit = T, plotreg = T) {
  
  if(class(modelList) != "list") modelList = list(modelList)
  
  # Separate variables
  vars = unlist(strsplit(deparse(.formula), "~"))
  
  y = gsub(" ", "", vars[1])
  
  x = gsub(" ", "", vars[2])
  
  # Extract model from modelList regressing the y variable
  
  y.model = modelList[[which(sapply(modelList, function(i) all.vars(formula(i))[1]) == y)]]
  
  if(is.null(y.model)) 
    
    stop("Check spelling of correlated variables - must match exactly response in model formula!")
  
  if(all(strsplit(deparse(formula(y.model)[[3]]), ".\\+.")[[1]] %in% x)) 
    
    stop("Y is a direct function of X, no partial residuals obtainable")

  # Get model formula 
  rhs = formula(drop.terms(y.model$terms, grep(x, attr(y.model$terms, "term.labels")), keep.response = TRUE))

  random.formula = get.random.formula(y.model, rhs, modelList, drop = x)
  
  # Get model control
  control = get.model.control(y.model, model.control)

  # Update model
  y.nox.model = suppressWarnings(if(is.null(random.formula)) 
    
    update(y.model, fixed = formula(rhs)) else
      
      if(any(class(y.model) %in% c("lme", "glmmPQL"))) 
        
        update(y.model, fixed = formula(rhs), random = formula(random.formula), control = control) else
          
          update(y.model, formula = formula(y, " ~ ", rhs, " + ", random.formula, sep = ""), control = control) )
  
  # Replace x variable as response in y model
  x.noy.model = suppressWarnings(if(any(class(y.model) %in% c("lme", "glmmPQL"))) 
    
    update(y.model, fixed = reformulate(deparse(formula(y.nox.model)[[3]]), response = x), random = as.formula(random.formula), control = control) else
      
      update(y.model, formula = reformulate(deparse(formula(y.nox.model)[[3]]), response = x.sub), control = control)
    
  )
  
  # Extract residuals from models and bind into data.frame
  y1 = data.frame(.id = 1:length(resid(y.nox.model)), resid(y.nox.model))
  
  x1 = data.frame(.id = 1:length(resid(x.noy.model)), resid(x.noy.model))
  
  resids.data = merge(y1, x1, by = ".id")[, -1]
  
  colnames(resids.data) = gsub(" ", "", c(y, x))
    
  # Plot results and regression line
  if(plotit == TRUE)
    
    plot(resids.data[, 1] ~ resids.data[ ,2], 
         xlab = ifelse(length(attr(y.nox.model$terms, "term.labels")) <= 3,
                       ifelse(length(attr(y.nox.model$terms, "term.labels")[-1]) > 1, paste(y, paste(attr(y.nox.model$terms, "term.labels")[-1], collapse=" + "), sep = " | "), paste(y)),
                       paste(y, "| others") ),
         ylab = ifelse(length(attr(x.noy.model$terms, "term.labels")) <= 3,
                       ifelse(length(attr(x.noy.model$terms, "term.labels")[-1]) > 1, paste(x, paste(attr(x.noy.model$terms, "term.labels")[-1], collapse=" + "), sep = " | "), paste(x)),
                       paste(x, "| others") )) 
    
  if(plotit == TRUE & plotreg == TRUE) abline(lm(resids.data[ ,1] ~ resids.data[ ,2]), col = "red", lwd = 2)
  
  return(resids.data)
    
}