partial.resid = function(.formula = y ~ x, modelList, model.control = NULL, plotit = TRUE, plotreg = TRUE) {
  
  if(any(class(modelList) != "list")) modelList = list(modelList)
  
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
  rhs = formula(drop.terms(terms(y.model), grep(x, attr(terms(y.model), "term.labels")), keep.response = TRUE))

  random.formula = get.random.formula(y.model, rhs, modelList, drop = x)
  
  # Get model control
  control = get.model.control(y.model, model.control)

  # Update model
  y.nox.model = suppressWarnings(if(is.null(random.formula)) 
    
    update(y.model, formula(rhs), control = control) else
      
      if(any(class(y.model) %in% c("lme", "glmmPQL"))) 
        
        update(y.model, fixed = formula(rhs), random = formula(random.formula), control = control) else
          
          update(y.model, paste(deparse(rhs), " + ", random.formula, collapse = ""), control = control) 
    
    )
  
  # Replace x variable as response in y model
  x.noy.model = suppressWarnings(if(is.null(random.formula)) {
    
    if(any(class(y.model) %in% c("glm"))) {
    
      # Try to update model
      mod = try(suppressWarnings(suppressMessages(
        update(y.model, 
               reformulate(deparse(formula(y.nox.model)[[3]]), response = x), 
               control = control)
        ) ), silent = TRUE)
      
      if(class(mod) == "try-error")
        
        update(y.model, 
               reformulate(deparse(formula(y.nox.model)[[3]]), response = x), 
               family = gaussian(), 
               na.action = na.omit,
               control = control) else mod 
      
    } else 
      
      update(y.model, 
             reformulate(deparse(formula(y.nox.model)[[3]]), response = x), 
             control = control)
    
  } else
      
      if(any(class(y.model) %in% c("lme", "glmmPQL"))) 
        
        update(y.model, fixed = reformulate(deparse(formula(y.nox.model)[[3]]), response = x), random = formula(random.formula), control = control) else
          
          update(y.model, reformulate(deparse(formula(y.nox.model)[[3]]), response = x), control = control) 
  
  )

  # Extract residuals from models
  y1 = data.frame(.id = names(resid(y.nox.model)), resid(y.nox.model))
  
  x1 = data.frame(.id = names(resid(x.noy.model)), resid(x.noy.model))
  
  # Merge residuals and store in a data.frame
  resids.data = merge(y1, x1, by = ".id", all = TRUE)[, -1]
  
  colnames(resids.data) = gsub(" ", "", c(y, x))
    
  # Plot results and regression line
  if(plotit == TRUE)
    
    plot(resids.data[, 1] ~ resids.data[ ,2], 
         ylab = ifelse(length(attr(terms(y.nox.model), "term.labels")) <= 3,
                       ifelse(length(attr(terms(y.nox.model), "term.labels")[-1]) > 1, paste(y, paste(attr(terms(y.nox.model), "term.labels")[-1], collapse=" + "), sep = " | "), paste(y)),
                       paste(y, "| others") ),
         xlab = ifelse(length(attr(terms(x.noy.model), "term.labels")) <= 3,
                       ifelse(length(attr(terms(x.noy.model), "term.labels")[-1]) > 1, paste(x, paste(attr(terms(x.noy.model), "term.labels")[-1], collapse=" + "), sep = " | "), paste(x)),
                       paste(x, "| others") )) 
    
  if(plotit == TRUE & plotreg == TRUE) abline(lm(resids.data[ ,1] ~ resids.data[ ,2]), col = "red", lwd = 2)
  
  return(resids.data)
    
}