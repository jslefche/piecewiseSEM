get.condind.ppa = function(modelList, data, corr.errors = NULL, add.vars = NULL, 
                            basis.set = NULL, disp.conditional = FALSE,
                            .progressBar = TRUE) {
  
  if(is.null(basis.set)) basis.set = get.basis.set(modelList, corr.errors, add.vars) else basis.set = basis.set
  
  if(.progressBar == T) pb = txtProgressBar(min = 0, max = length(basis.set), style = 3) else pb = NULL
  
  pvalues.df = do.call(rbind, lapply(seq_along(basis.set), function(i) {
    
    basis.mod = modelList[[match(basis.set[[i]][2], unlist(lapply(modelList, function(j) as.character(formula(j)[2]))))]]
      
      #### Need to fix getting random effects structure    
      fixed.formula = rbind(
        if(length(basis.set[[i]])<3) paste(basis.set[[i]][2], "~", paste(basis.set[[i]][1])),
        if (length(basis.set[[i]])>2) paste(basis.set[[i]][2], "~", paste(basis.set[[i]][c(1, 3:length(basis.set[[i]]))], collapse = "+"))
      )
    
    
    random.formula = NULL
    
    modelList.random.slopes = NULL 
   
    random.formula = NULL
  
    if(is.null(random.formula)) basis.mod = update(basis.mod, fixed.formula, data = data) 
 
    if (!grepl(":|\\*", basis.set[[i]][1])) 
      if(all(class(basis.mod) %in% c("pgls"))) row.num = which(basis.set[[i]][1] ==  attr(coef(basis.mod), "names")) else
        row.num = which(basis.set[[i]][1] == attr(terms(basis.mod), "term.labels")) + 1  else 
          if(all(class(basis.mod) %in% c("pgls"))) row.num = which(grepl(":|\\*",  attr(coef(basis.mod), "names"))) else
            row.num = which(grepl(":|\\*", attr(terms(basis.mod), "term.labels"))) + 1 
    
    ret = if(all(class(basis.mod) %in% c("lm", "glm", "negbin", "glmerMod", "merModLmerTest","pgls"))) 
      as.data.frame(t(unname(summary(basis.mod)$coefficients[row.num, ]))) else
        as.data.frame(t(unname(summary(basis.mod)$tTable[row.num, ])))
    
    if(length(ret) != 5) ret = cbind(ret[1:2], NA, ret[3:4])
    
    names(ret)=c("estimate","std.error","DF","crit.value","p.value")
    
    
    if(disp.conditional == TRUE) 
      ret = cbind(missing.path = paste(basis.set[[i]][2], "<-", paste(basis.set[[i]][1], collapse = "+")), 
                  conditional.on = paste(basis.set[[i]][3:length(basis.set[[i]])], collapse = ","),
                  ret) else
                    ret = cbind(missing.path = paste(basis.set[[i]][2], "<-", paste(basis.set[[i]][1], collapse = "+")), 
                                ret)
    
    if(.progressBar == TRUE) setTxtProgressBar(pb, i)
    
    return(ret)
    
  } ) )
  
  if(!is.null(pb)) close(pb)  
  
  return(pvalues.df)
  
}