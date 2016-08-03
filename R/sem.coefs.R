sem.coefs = function(modelList, data = NULL, standardize = "none", corr.errors = NULL, intercept = FALSE) {
  
  if(any(class(modelList) != "list")) modelList = list(modelList)
  
  names(modelList) = NULL

  if(!standardize %in% c("none", "scale", "range")) stop("'standardize' must equal 'none', 'scale', or 'range'")
  
  # Scale variables, if indicated
  if(standardize != "none") newdata = get.scaled.data(modelList, data, standardize)
  
  # Return coefficients
  ret = do.call(rbind, lapply(modelList, function(i) {
    
    if(standardize != "none") i = get.scaled.model(i, newdata, modelList)

    if(intercept == TRUE) irow = TRUE else irow = -1
    
    # Extract coefficients and return in a data.frame
    if(any(class(i) %in% c("lm", "glm", "pgls", "negbin", "glmerMod", "glmmadmb"))) {
      
      tab = summary(i)$coefficients
      
      data.frame(response = Reduce(paste, deparse(formula(i)[[2]])),
                 predictor = rownames(tab)[irow],
                 estimate = tab[irow, 1],
                 std.error = tab[irow, 2],
                 p.value = tab[irow, 4]
                 )
      
    } else if(any(class(i) %in% c("gls"))) {
      
      tab = summary(i)$tTable
      
      data.frame(response = Reduce(paste, deparse(formula(i)[[2]])),
                 predictor = rownames(tab)[irow],
                 estimate = tab[irow, 1],
                 std.error = tab[irow, 2],
                 p.value = tab[irow, 4]
                 )
      
    } else if(any(class(i) %in% c("lme", "glmmPQL"))) {
      
      tab = summary(i)$tTable
      
      data.frame(response = Reduce(paste, deparse(formula(i)[[2]])),
                 predictor = rownames(tab)[irow],
                 estimate = tab[irow, 1],
                 std.error = tab[irow, 2],
                 p.value = tab[irow, 5]
                 )
      
    } else if(any(class(i) %in% c("lmerMod", "merModLmerTest"))) {
  
      tab = suppressMessages(summary(i)$coefficients)
          
      # Loop over variables and return P-values
      kr.p = sapply(names(fixef(i))[irow], function(x) {
        
        i.reduced = update(i, as.formula(paste("~ . -", x)))
      
        KRmodcomp(i, i.reduced)$test$p.value[1]
        
      } )
      
      # KRSumFun <- function(object, objectDrop, ...) {
      #   krnames <- c("ndf","ddf","Fstat","p.value","F.scaling")
      #   r <- if (missing(objectDrop)) {
      #     setNames(rep(NA,length(krnames)),krnames)
      #   } else {
      #     krtest <- KRmodcomp(object,objectDrop)
      #     unlist(krtest$stats[krnames])
      #   }
      #   attr(r,"method") <- c("Kenward-Roger via pbkrtest package")
      #   r
      # }
      # 
      # kr.p = drop1(i, test = "user", sumFun = KRSumFun)

      data.frame(response = Reduce(paste, deparse(formula(i)[[2]])),
                 predictor = rownames(tab)[irow],
                 estimate = tab[irow, 1],
                 std.error = tab[irow, 2],
                 p.value = kr.p
                 ) 
      } 
    
    } ) )
  
  # Do significance tests for correlated errors
  if(!is.null(corr.errors)) 
    
    ret = rbind(ret, do.call(rbind, lapply(corr.errors, function(j) {
      
      # Pull out correlated variables
      corr.vars = gsub(" ", "", unlist(strsplit(j, "~~")))
      
      # Final model with response
      corr.mod = modelList[[match(corr.vars[1], sapply(modelList, function(k) paste(formula(k)[2])))]]
      
      if(!is.null(corr.mod)) {
        
        # Update model to include correlated error
        if(any(class(corr.mod) %in% c("lme", "glmmPQL")))
          
          corr.mod = update(corr.mod, fixed = formula(paste0(". ~ . + ", corr.vars[2]))) else
            
            corr.mod = update(corr.mod, formula(paste0(". ~ . + ", corr.vars[2])))
        
        # Get partial residuals
        corr.mod.resids = partial.resid(formula(paste0(corr.vars, collapse = " ~ ")), corr.mod, data, plotit = FALSE)
        
        # Perform significance test and return in a data.frame
        data.frame(
          response = paste("~~", corr.vars[1]),
          predictor = paste("~~", corr.vars[2]),
          estimate = cor(corr.mod.resids[, 1], corr.mod.resids[, 2], use = "complete.obs"),
          std.error = NA,
          p.value =  1 - pt(
            (cor(corr.mod.resids[, 1], corr.mod.resids[, 2], use = "complete.obs") * sqrt(nrow(data) - 2))/
                 (sqrt(1 - cor(corr.mod.resids[, 1], corr.mod.resids[, 2], use = "complete.obs")^2)), nrow(data)-2)
          )
        
      } else {
        
        data.frame(
          response = paste("~~", corr.vars[1]),
          predictor = paste("~~", corr.vars[2]),
          estimate = cor(data[, corr.vars[1]], 
                         data[, corr.vars[2]], 
                         use = "complete.obs"),
          std.error = NA,
          p.value = 1 - 
            pt((cor(data[, corr.vars[1]], data[, corr.vars[2]], use = "complete.obs") * sqrt(nrow(data) - 2))/
                 (sqrt(1 - cor(data[, corr.vars[1]], data[, corr.vars[2]], use = "complete.obs")^2)), nrow(data)-2)
          )
        
      }
      
    } ) ) )
  
  # Order by response and p-value
  ret = ret[with(ret, order(response, p.value)),]
  
  # Round all numeric values
  # ret[, sapply(ret, is.numeric)] = apply(ret[, sapply(ret, is.numeric)], 2, round, 4)
  
  # Round only P-values
  ret$p.value = round(ret$p.value, 4)
  
  # Assign significance indicators
  sig = sapply(ret$p.value, function(y) {
    
    ifelse(y > 0.01 & y < 0.05, "*", 
           ifelse(y > 0.001 & y <= 0.01, "**",
                  ifelse(y <= 0.001, "***", "")
           )
    )
    
  } )
  
  ret = cbind(ret, sig)
  
  colnames(ret)[ncol(ret)] = ""
  
  # Remove rownames
  rownames(ret) = NULL
  
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