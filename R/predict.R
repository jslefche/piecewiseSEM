#' @title Generate predictions from a fit piecewise SEM
#'
#' @param object An object of class [psem]
#' @param newdata Optional. A data frame with new values
#' of exogenous variables
#' @param se.fit Do you want to incorporate coefficient 
#' standard errors? Currently defaults to FALSE as it has
#' not been implemented yet.
#' @param interval What type of interval do we want to 
#' report - fit, prediction, or none. Currently not 
#' implemented, so "none" is defaulted. Arugment here for
#' later development.
#' @param ... Other arguments to `predict` functions used by
#' pieces of the `psem` object. 
#'
#' @return A data.frame of predicted values
#' @export
#'
#' @examples
#' data(keeley)
#' fit <- psem(
#'   lm(abiotic ~ distance, data=keeley),
#'   lm(rich ~ abiotic + hetero, data=keeley),
#'   lm(hetero ~ distance, data=keeley),
#'   data = keeley)
#' 
#' fit_vals <- predict(fit)
#' 
#' new_predictions <- 
#'   predict(fit, 
#'           newdata=data.frame(distance=c(30, 50)))
#'           
#' head(fit_vals)
#' head(new_predictions)

predict.psem <- function(object, newdata=NULL, 
                                 se.fit=FALSE,
                                 interval = "none",
                                 ...){
  
  #get sorted psem object to work with
  object <- getSortedPsem(object)
  object <- removeData(object, formulas = 1)
  
  #Now, get formulae
  formulaList <- listFormula(object)
  lhs <- getLHS(formulaList)
  
  #If there is no desire for intervals
  if(se.fit==FALSE && interval=="none"){
    
    #If there is no new data, just get predict from each function  
    if(is.null(newdata)){
      ret <- predict_piecewiseSEM_nodata(object, ...)
    }else{
      
      #otherwise, some simple looping will suffice
      ret <- predict_piecewiseSEM_newdata(object, newdata, lhs, ...)
      
    }
    
  }else{
    stop("Error intervals not yet implemented.")
    #If there is no new data, make new data with
    #just the exogenous variables
    #rhs <- getRHS(formulaList)
    #exo <- rhs %not_in% lhs
    
    #Now, get predictions from propagated simulations
  }
  
  # This this was fixed with the psem sorting - but keeping
  # in case not
  # #fix column names
  # if(length(lhs)==ncol(ret)) {
  #   names(ret) <- lhs
  # }else{
  #   lhs_new <- do.call(c,lapply(lhs, rep, ncol(ret)/length(lhs), simplify=TRUE))
  #   names(ret) <- paste(lhs_new, names(ret), sep="_")
  # }
  
  return(ret)
}


#' @title Get fit values from a fitted piecewise SEM
#'
#' @param object A  [psem] object
#'
#' @return A data.frame of fitted values
#' @export
#'
#' @examples
#' data(keeley)
#' fit <- psem(
#'   lm(abiotic ~ distance, data=keeley),
#'   lm(rich ~ abiotic + hetero, data=keeley),
#'   lm(hetero ~ distance, data=keeley),
#'   data = keeley)
#' 
#' fit_vals <- fitted(fit)
#'           
#' head(fit_vals)

fitted.psem <- function(object){
  predict.psem(object)
}

#' Extract fit values from models making up psem
#' 
#' @keywords internal
#' 
predict_piecewiseSEM_nodata <- function(modList, ...){
  ret <- lapply(modList, predict, ...)
  ret <- as.data.frame(do.call(cbind, ret))
  ret
}

#' Generate fit values newdata with psem model
#' 
#' @keywords internal
#' 
predict_piecewiseSEM_newdata <- function(modList, newdata, lhs_sorted, ...){
  #If there is new data...
  modData <- newdata
  for(j in 1:length(modList)){
    pred <-  as.data.frame(predict(modList[[j]], newdata=modData, ...))
    modData <- cbind(modData, pred[,1])
    names(modData)[j+1] <- lhs_sorted[j]
  }
  ret <- modData[,-1]
}

