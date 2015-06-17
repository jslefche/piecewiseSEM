sem.coefs = function(modelList, data, standardize = "none", corr.errors = NULL) {
  
  if(any(class(modelList) != "list")) modelList = list(modelList)
  
  names(modelList) = NULL

  if(!standardize %in% c("none", "scale", "range")) stop("'standardize' must equal 'none', 'scale', or 'range.'")
  
  # Scale variables, if indicated
  if(standardize != "none") {
    
    # Get variables to scale, ignoring variables that are modeled to non-normal distributions
    vars.to.scale = unlist(lapply(modelList, function(i) {
     
      err = try(family(i), TRUE)
      
      if(grepl("Error", err[1]) | grepl("gaussian", err[1]))
        
        rownames(attr(terms(i), "factors")) else {
          
          print(
            paste("Reponse '", formula(i)[2], "' is not modeled to a gaussian distribution: keeping response on original scale")
          )
          
          NULL }

      }
    
    ) )
    
    # Remove variables that are factors
    vars.to.scale = vars.to.scale[!vars.to.scale %in% colnames(data)[sapply(data, is.factor)]]
    
    # Remove duplicated variables
    vars.to.scale = vars.to.scale[!duplicated(vars.to.scale)]
    
    # Scale those variables by mean and SD, or by range
    data[, vars.to.scale] = apply(data[, vars.to.scale], 2, function(x) 
      
      if(standardize == "scale") scale(x) else
        
        if(standardize == "range") (x-min(x, na.rm = T)) / diff(range(x, na.rm = T))
      
      )
    
    }
  
  # Return coefficients
  ret = do.call(rbind, lapply(modelList, function(i) {
    
    if(standardize != "none") i = update(i, data = data)
    
    # Extract coefficients and return in a data.frame
    if(any(class(i) %in% c("lm", "glm", "pgls", "negbin", "glmerMod"))) {
      
      tab = summary(i)$coefficients
      
      data.frame(response = Reduce(paste, deparse(formula(i)[[2]])),
                 predictor = rownames(tab)[-1],
                 estimate = tab[-1, 1],
                 std.error = tab[-1, 2],
                 p.value = tab[-1, 4], 
                 row.names = NULL)
      
    } else if(any(class(i) %in% c("gls"))) {
      
      tab = summary(i)$tTable
      
      data.frame(response = Reduce(paste, deparse(formula(i)[[2]])),
                 predictor = rownames(tab)[-1],
                 estimate = tab[-1, 1],
                 std.error = tab[-1, 2],
                 p.value = tab[-1, 4], 
                 row.names = NULL)
      
    } else if(any(class(i) %in% c("lme", "glmmPQL"))) {
      
      tab = summary(i)$tTable
      
      data.frame(response = Reduce(paste, deparse(formula(i)[[2]])),
                 predictor = rownames(tab)[-1],
                 estimate = tab[-1, 1],
                 std.error = tab[-1, 2],
                 p.value = tab[-1, 5], 
                 row.names = NULL)
      
    } else if(any(class(i) %in% c("lmerMod", "merModLmerTest"))) {
      
      tab = summary(as(i, "merModLmerTest"))$coefficients
      
      data.frame(response = Reduce(paste, deparse(formula(i)[[2]])),
                 predictor = rownames(tab)[-1],
                 estimate = tab[-1, 1],
                 std.error = tab[-1, 2],
                 p.value = tab[-1, 5], 
                 row.names = NULL) 
      
      } 
    
    } ) )
  
  # Do significance tests for correlated errors
  if(!is.null(corr.errors)) 
    
    ret = rbind(ret, do.call(rbind, lapply(corr.errors, function(j) {
      
      # Pull out correlated variables
      corr.vars = gsub(" ", "", unlist(strsplit(j, "~~")))
      
      # Perform significance test and return in a data.frame
      data.frame(
        response = paste("~~", corr.vars[1]),
        predictor = paste("~~", corr.vars[2]),
        estimate = cor(data[, corr.vars[1]], 
                       data[, corr.vars[2]], 
                       use = "complete.obs"),
        std.error = NA,
        p.value = 1 - 
          pt((cor(data[, corr.vars[1]], data[, corr.vars[2]], use = "complete.obs") * sqrt(nrow(data) - 2))/
               (sqrt(1 - cor(data[, corr.vars[1]], data[, corr.vars[2]], use = "complete.obs")^2)), nrow(data)-2),
        row.names = NULL
        )
      
    } ) ) )
  
  # Order by p-value
  ret = ret[order(ret$p.value), ]
  
  return(ret)

}