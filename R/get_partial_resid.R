get.partial.resid = function(.formula = y ~ x, modelList, model.control = NULL, plotit= T, regr = T) {
  
  if(class(modelList) != "list") modelList = list(modelList)
  
  vars = unlist(strsplit(deparse(.formula), "~"))
  y = gsub(" ", "", vars[1])
  x = gsub(" ", "", vars[2])
  
  y.model = modelList[[match(y, unlist(lapply(modelList, function(i) gsub(".*\\((.*?)\\)+.*", "\\1", formula(i)[2]))))]]
  
  if(is.null(y.model)) stop("Check spelling of correlated variables - must match exactly response in model formula!")
  
  if(all(x %in% strsplit(deparse(formula(y.model)[[3]]), "\\+")[[1]])) 
    stop("Y is a direct function of X, no partial residuals obtainable")
  
  if(length(model.control)>10) model.control = list(model.control)
  
  if(is.null(model.control)) {
    if(class(y.model) %in% c("lme", "glmmPQL")) control = lmeControl() else 
      if(class(y.model) %in% c("lmerMod", "merModLmerTest")) control = lmerControl() else
        if(class(y.model) %in% c("glmerMod")) control = glmerControl() else
          control = NULL    
  } else {
    if(!is.null(model.control) & class(y.model) %in% c("lme", "glmmPQL"))
      control = model.control[[which(sapply(lapply(model.control, function(x) attr(x, "class")), is.null))]] else
        if(!is.null(model.control) & class(y.model) %in% c("lmerMod", "merModLmerTest"))
          control = model.control[[which(sapply(lapply(model.control, function(x) attr(x, "class")), function(x) any(x %in% "lmerControl")))]] else
            if(!is.null(model.control) & class(y.model) %in% c("glmerMod"))
              control = model.control[[which(sapply(lapply(model.control, function(x) attr(x, "class")), function(x) any(x %in% "glmerControl")))]] else
                control = NULL }
  
  if(class(y.model) %in% c("lme", "glmmPQL")) {
    if(!grepl("\\+", deparse(y.model$call$random))) random = y.model$call$random else {
      random = drop.terms(terms(as.formula(gsub("\\|", "+BERLIOZ+", deparse(y.model$call$random))), special = c("+BERLIOZ+","\\/")),
                          which(all.vars(y.model$call$random) == x))
      random = as.formula(gsub("\\+.BERLIOZ.\\+", "\\|", deparse(random)))
      if(grepl("\\:", deparse(random))) random = as.formula(gsub("..[^\\+]+:", "\\/", deparse(random))) }
    y.nox.model = update(y.model, fixed = drop.terms(y.model$terms, grep(x, attr(y.model$terms, "term.labels")), keep.response = TRUE), random = as.formula(random), control = control) 
  } else {
    y.nox.model = suppressWarnings(
      update(y.model, formula = drop.terms(y.model$terms, grep(x, attr(y.model$terms, "term.labels")), keep.response = TRUE), control = control) )
  }
  
  x.sub = attr(y.model$terms, "term.labels")[grep(x, attr(y.model$terms, "term.labels"))][1]
  
  if(all(class(y.nox.model) %in% c("lme", "glmmPQL"))) 
    x.noy.model = update(y.model, fixed = reformulate(deparse(formula(y.nox.model)[[3]]), response = x.sub), random = as.formula(random), control = control) else 
      x.noy.model = suppressWarnings( 
        update(y.model, formula = reformulate(deparse(formula(y.nox.model)[[3]]), response = x.sub), control = control) )
    
  y1 = data.frame(.id = 1:length(resid(y.nox.model)), resid(y.nox.model))
  x1 = data.frame(.id = 1:length(resid(x.noy.model)), resid(x.noy.model))
  resids.data = merge(y1, x1, by = ".id")[, -1]
  
  colnames(resids.data) = gsub(" ", "", c(y, x))
    
  if(plotit == TRUE) plot(resids.data[, 1] ~ resids.data[ ,2], 
                          xlab = ifelse(length(attr(y.nox.model$terms, "term.labels")) <= 3,
                                        ifelse(length(attr(y.nox.model$terms, "term.labels")[-1]) > 1, paste(y, paste(attr(y.nox.model$terms, "term.labels")[-1], collapse=" + "), sep = " | "), paste(y)),
                                        paste(y, "| others") ),
                          ylab = ifelse(length(attr(x.noy.model$terms, "term.labels")) <= 3,
                                        ifelse(length(attr(x.noy.model$terms, "term.labels")[-1]) > 1, paste(x, paste(attr(x.noy.model$terms, "term.labels")[-1], collapse=" + "), sep = " | "), paste(x)),
                                        paste(x, "| others") )) 
    
  if(plotit == TRUE & regr == TRUE) abline(lm(resids.data[ ,1] ~ resids.data[ ,2]), col = "red", lwd = 2)
  
  return(resids.data)
    
}