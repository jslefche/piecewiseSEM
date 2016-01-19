sem.coefs = function(modelList, data, standardize = "none", corr.errors = NULL) {
  
  if(any(class(modelList) != "list")) modelList = list(modelList)
  
  names(modelList) = NULL

  if(!standardize %in% c("none", "scale", "range")) stop("'standardize' must equal 'none', 'scale', or 'range'")
  
  # Scale variables, if indicated
  if(standardize != "none") newdata = scale.data(modelList, data, standardize)
  
  # Return coefficients
  ret = do.call(rbind, lapply(modelList, function(i) {
    
    if(standardize != "none") i = update(i, data = newdata)
    
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
  
  # Order by response and p-value
  ret = ret[with(ret, order(response, p.value)),]
  
  # Round p-values
  ret$p.value = round(ret$p.value, 4)
  
#   # If standardize != "none" and interactions present, set SEs and P-values to NA
#   if(standardize != "none" & any(sapply(modelList, function(x) any(grepl("\\:|\\*", formula(x)))))) {
#     
#     # Return warning
#     print("It is not correct to interpret significance of standardized variables involved in interactions!")
#     print("Refer to unstandardized P-values to assess significance.") 
#   
#     # Remove SEs and P=values for rows with interactions
#     ret[grepl("\\:|\\*", ret$predictor), 4:5] = NA
#     
#   }
    
  return(ret)

}