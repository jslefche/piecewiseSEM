get.scaled.model = function(model, newdata) {
  
  if(any(class(model) %in% c("lmerMod", "merModLmerTest", "glmerMod"))) {
    
    # Get random effects
    rand.effs = gsub(" ", "", sapply(findbars(formula(model)), function(x) gsub(".*\\|(.*)", "\\1", deparse(x))))
    
    # Get fixed effects
    fixed.effs = all.vars(formula(model))[!all.vars(formula(model)) %in% rand.effs]
    
    # Get fixed formula stripped of transformations
    fixed.form = paste0(fixed.effs[1], " ~ ", paste0(fixed.effs[-1], collapse = " + "))
    
    # Bind back in random structure
    random.form = get.random.formula(model, rhs = paste0(fixed.effs[-1], collapse = " + "), modelList)
    
    # Get updated formula
    new.form = paste0(fixed.form, " + ", random.form)
    
  } else {
    
    new.form = paste0(all.vars(formula(model))[1], " ~ ", paste0(all.vars(formula(model))[-1], collapse = " + "))
    
  }
  
  # Update model
  if(any(class(model) %in% c("lme", "glmmPQL")))
    
    model = update(model, fixed = as.formula(new.form), data = newdata) else
      
      model = update(model, formula = new.form, data = newdata)
    
  # Return model
  return(model)
  
}