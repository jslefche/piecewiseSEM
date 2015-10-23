partial.resid = function(
  
  .formula = y ~ x, modelList, data, model.control = NULL, return.data.frame = TRUE, plotit = TRUE, plotreg = TRUE, plotCI = TRUE
  
  ) {
  
  if(any(class(modelList) != "list")) modelList = list(modelList)
  
  # Separate variables
  vars = unlist(strsplit(deparse(.formula), "~"))
  
  y = gsub(" ", "", vars[1])
  
  x = gsub(" ", "", vars[2])
  
  # Extract model from modelList regressing the y variable
  y.vars = suppressWarnings(sapply(modelList, function(i) all.vars(formula(i))[1]) == y)
  
  if(!any(y.vars)) 
    
    stop("Check spelling of correlated variables - must match exactly response in model formula!") else
      
      y.model = modelList[[which(y.vars)]]
  
  if(all(strsplit(deparse(formula(y.model)[[3]]), ".\\+.")[[1]] %in% x)) 
    
    stop("Y is a direct function of X, no partial residuals obtainable")

  # Get model formula 
  rhs = formula(drop.terms(terms(y.model), grep(gsub("\\*", "\\:", x), attr(terms(y.model), "term.labels")), keep.response = TRUE))

  random.formula = get.random.formula(y.model, rhs, modelList, dropterms = x)
  
  # Get model control
  control = get.model.control(y.model, model.control)

  # Update model
  y.nox.model = suppressWarnings(if(is.null(random.formula)) 
    
    update(y.model, formula(rhs), control = control) else
      
      if(any(class(y.model) %in% c("lme", "glmmPQL"))) 
        
        update(y.model, fixed = formula(rhs), random = random.formula, control = control) else
          
          update(y.model, formula(paste(rhs, " + ", random.formula, collapse = "")), control = control) 
    
    )
  
  # If x is an interaction, calculate and replace in dataset
  if(grepl("\\:|\\*", x)) { 
    
    data[, gsub("\\:|\\*", "_", x)] = apply(data[, strsplit(x, "\\:|\\*")[[1]]], 1, prod)
    
    x = gsub("\\:|\\*", "_", x)
    
  }
  
  # Replace x variable as response in y model
  x.noy.model = suppressWarnings(if(is.null(random.formula)) {
    
    if(any(class(y.model) %in% c("glm"))) {
    
      # Try to update model
      mod = try(suppressWarnings(suppressMessages(
        update(y.model, 
               reformulate(deparse(formula(y.nox.model)[[3]]), response = x), 
               control = control, 
               data = data)
        ) ), silent = TRUE)
      
      if(class(mod) == "try-error")
        
        update(y.model, 
               reformulate(deparse(formula(y.nox.model)[[3]]), response = x), 
               family = gaussian(), 
               na.action = na.omit,
               control = control,
               data = data) else 
                 
                 mod 
      
    } else 
      
      update(y.model, 
             reformulate(deparse(formula(y.nox.model)[[3]]), response = x), 
             control = control,
             data = data)
    
  } else
      
      if(any(class(y.model) %in% c("lme", "glmmPQL"))) 
        
        update(y.model, fixed = reformulate(deparse(formula(y.nox.model)[[3]]), response = x), 
               random = random.formula, 
               control = control,
               data = data) else
          
          update(y.model, 
                 reformulate(deparse(formula(y.nox.model)[[3]]), response = x), 
                 control = control,
                 data = data) 
  
  )

  # Extract residuals from models
  if(any(class(y.nox.model) %in% c("lme", "glmmPQL"))) {
    
    # Get innermost level of grouping
    Q = length(summary(y.nox.model)$modelStruct$reStruct)
    
    y.resids = resid(y.nox.model, level = 0:Q)
    
    y.resids = y.resids[, 0, drop = FALSE]
    
    } else y.resids = resid(y.nox.model)
  
  if(any(class(x.noy.model) %in% c("lme", "glmmPQL"))) {
    
    # Get innermost level of grouping
    Q = length(summary(x.noy.model)$modelStruct$reStruct)
    
    x.resids = resid(x.noy.model, level = 0:Q) 
    
    x.resids = x.resids[, 0, drop = FALSE]
    
    } else x.resids = resid(x.noy.model)
  
  # Bind together in data.frame
  y1 = data.frame(.id = names(y.resids), y.resids)
  
  x1 = data.frame(.id = names(x.resids), x.resids)
  
  # Merge residuals and store in a data.frame
  resids.data = merge(y1, x1, by = ".id", all = TRUE)[, -1]
  
  colnames(resids.data) = gsub(" ", "", c(y, x))
    
  # Plot results and regression line
  if(plotit == TRUE)
    
    plot(resids.data[, 1] ~ resids.data[ ,2], 
         ylab = ifelse(nchar(paste(attr(terms(y.nox.model), "term.labels"), y,  collapse = " + ")) <= 20,
                       paste(y, paste(attr(terms(y.nox.model), "term.labels"), collapse=" + "), sep = " | "), 
                       paste(y, "| others") ),
         xlab = ifelse(nchar(paste(attr(terms(x.noy.model), "term.labels"), x,  collapse = " + ")) <= 20,
                       paste(gsub("_", "\\*", x), paste(attr(terms(x.noy.model), "term.labels"), collapse=" + "), sep = " | "),
                       paste(gsub("_", "\\*", x), "| others") )
    )

  if(plotit == TRUE & plotreg == TRUE | plotCI == TRUE) {
    
    # Get names of resids.data
    resids.y = names(resids.data)[1]
    resids.x = names(resids.data)[2]
    
    # Regression y ~ x
    new.mod = lm(formula(paste(resids.y, " ~ ", resids.x, sep = "")), resids.data)
    
    # Create new data.frame
    newdata = data.frame(
      seq(min(resids.data[, 2], na.rm = TRUE) - min(resids.data[, 2], na.rm = TRUE) * 0.1, 
          max(resids.data[, 2], na.rm = TRUE) + max(resids.data[, 2], na.rm = TRUE) * 0.1, 
          length.out = nrow(resids.data) * 2)
    )
    
    colnames(newdata) = resids.x
    
    # Generate predictions
    pred = predict(new.mod, newdata, interval = "confidence", level = 0.95)
    
    # Plot fitted curve
    if(plotit == TRUE) abline(new.mod, col = "red", lwd = 2)
    
    # Plot confidence intervals
    if(plotCI == TRUE) { 
      
      lines(newdata[, 1], pred[, 2], col = "red", lwd = 1.8, lty = 2)
      
      lines(newdata[, 1], pred[, 3], col = "red", lwd = 1.8, lty = 2)
    
    }
    
  }
  
  if(return.data.frame == TRUE) return(resids.data)
    
}