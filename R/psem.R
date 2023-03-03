#' Fitting piecewise structural equation models
#'
#' \code{psem} is used to unite a list of structural equations into a single
#' structural equation model.
#'
#' \code{psem} takes a list of structural equations, which can be model objects
#' of classes: \code{lm, glm, gls, pgls, Sarlm, lme, glmmPQL, lmerMod,
#' lmerModLmerTest, glmerMod, glmmTMB, gam}.
#'
#' It also takes objects of class \code{formula, formula.cerror}, corresponding
#' to additional variables to be included in the tests of directed separation
#' (\code{X ~ 1}) or correlated errors (\code{X1 \%~~\% X2}).
#'
#' The function optionally accepts data objects of classes: \code{matrix,
#' data.frame, SpatialPointsDataFrame, comparative.data}, or these are derived
#' internally from the structural equations.
#'
#' @param \dots A list of structural equations
#' 
#' @return Returns an object of class \code{psem}
#' 
#' @author Jon Lefcheck <LefcheckJ@@si.edu>
#' @seealso \code{\link{summary.psem}}, \code{\link{\%~~\%}}
#' 
#' @examples 
#' mod <- psem(
#' lm(rich ~ cover, data = keeley),
#' lm(cover ~ firesev, data = keeley),
#' lm(firesev ~ age, data = keeley),
#' data = keeley
#' )
#' 
#' summary(mod)
#' 
#' @export
#' 
psem <- function(...) {

  x <- list(...)

  x <- formatpsem(x)

  class(x) <- "psem"
  
  x
  
}

#' Format for psem
#' 
#' @keywords internal
#' 
formatpsem <- function(x) {

  evaluateClasses(x)
  
  stop_psem(x)
  
  # checkTransformations(x)
  
  idx <- which(sapply(x, function(y) any(class(y) %in% c("matrix", "data.frame", "SpatialPointsDataFrame", "comparative.data"))))

  if(sum(idx) == 0) idx <- which(names(x) == "data")

  if(sum(idx) > 0) {

    if(is.null(names(x))) names(x) <- 1:length(x)

    names(x)[idx] <- "data"

    x$data <- x$data

  } else {

    x$data <- GetData(x)

  }
  
  x <- checkData(x)

  if(any(is.na(names(x)))) {

    idx. <- which(is.na(names(x)))

    names(x)[idx.] <- idx.

  }
  
  # if(!any(sapply(x, class) %in% "gam")) {
  # 
  #   vars <- as.vector(unlist(sapply(removeData(x, formulas = 1), function(x) unname(all.vars_notrans(x)))))
  # 
  #   vars <- vars[!duplicated(vars) & !grepl("\\:", vars)]
  #   
  #   t_vars <- as.vector(unlist(sapply(removeData(x, formulas = 1), function(x) unname(all.vars_trans(x, smoothed = TRUE)))))
  #   
  #   t_vars <- t_vars[!duplicated(t_vars) & !grepl("\\:", t_vars)]
  #   
  #   if(length(vars) != length(t_vars)) stop("Some variables appear as alternately transformed and untransformed. Apply transformations across the entire model", call. = FALSE)
  #   
  # }
    
  if(all(class(x$data) == "comparative.data")) { 
    
    if(any(sapply(x$data$data, is.na))) warning("NAs detected in the dataset. Consider removing all rows with NAs to prevent fitting to different subsets of data", call. = FALSE) 
    
    } else {
      
      if(any(sapply(x$data, is.na))) warning("NAs detected in the dataset. Consider removing all rows with NAs to prevent fitting to different subsets of data", call. = FALSE)

    }

  formulaList <- listFormula(x, formulas = 1)

  if(any(duplicated(sapply(formulaList, function(y) all.vars.merMod(y)[1]))))

    stop("Duplicate responses detected in the model list. Collapse into single multiple regression!", call. = FALSE)

  x

}

#' Convert list to psem object
#' 
#' @param object any \code{R} object
#' @param Class the name of the class to which \code{object} should be coerced
#' 
#' @export
#' 
as.psem <- function(object, Class = "psem") { 
  
  object <- formatpsem(object)
  
  class(object) <- Class
  
  object
  
}

#' Check to see whether supplied data.frame matches model-extracted data
#' 
#' @keywords internal
#' 
checkData <- function(x) {
  
  if(!identical(x$data, GetData(removeData(x)))) x$data <- GetData(removeData(x))
  
  return(x)
  
}

#' Check to see whether variables exist as transformed and untransformed
#' 
#' @keywords internal
#' 
checkTransformations <- function(x) {
  
  
}

#' Evaluate model classes and stop if unsupported model class
#' 
#' @param x a list of structural equations or a model object
#' 
#' @export
#' 
evaluateClasses <- function(x) {

  classes <- unlist(sapply(x, class))

  classes <- classes[!duplicated(classes)]

  supported.classes <- c(
    "character",
    "matrix", "data.frame", "SpatialPointsDataFrame", "comparative.data",
    "data.table", "tbl_df", "tbl",
    "formula", "formula.cerror",
    "lm", "glm", "gls", "negbin",
    "lme", "glmmPQL",
    "lmerMod", "lmerModLmerTest", "glmerMod", "glmmTMB",
    "Sarlm",
    "pgls", "phylolm", "phyloglm",
    "gam"
  )

  if(!all(classes %in% supported.classes))

    stop(
      paste0(
        "Unsupported class in model list: ",
        paste0(classes[!classes %in% supported.classes], collapse = ", "),
        ". See 'help(piecewiseSEM)' for more details."),
      call. = FALSE
    )

}

#' Stop function for unsupported methods
#' 
#' @keywords internal
#' 
stop_psem <- function(x) {
  
  if(any(sapply(x, function(y) grepl("poly\\(.*\\)", formula(y))))) stop("Polynomials not supported", call. = F)
  
}

#' Print psem
#' 
#' @param x an object of class psem
#' @param ... further arguments passed to or from other methods
#' 
#' @method print psem
#' 
#' @export
#' 
print.psem <- function(x, ...) {

  nm <- deparse(substitute(x))
  
  formulas <- listFormula(x)

  formulas_print <- sapply(1:length(formulas), function(i) {

    if(inherits(formulas[[i]], "formula.cerror"))

      paste0("Correlated error: ", paste(formulas[[i]])) else

        paste0(class(x[[i]])[1], ": ", paste0(deparse(formulas[[i]]), collapse = ""))

  } )

  data_print <- if(!is.null(x$data)) head(x$data) else head(GetData(x))

  class_print <- paste0("class(", class(x), ")")

  cat("Structural Equations of", nm, ":\n")

  cat(paste(formulas_print, collapse = "\n"))

  cat("\n\nData:\n")

  print(data_print)

  cat(paste("...with ", dim(x$data)[1] - 6, " more rows"))

  cat("\n\n")

  print(class_print)

}

#' Update psem model object with additional values.
#' 
#' @param object a psem object
#' @param ... additional arguments to update
#' 
#' @method update psem
#' 
#' @examples 
#' mod <- psem(
#' lm(rich ~ cover, data = keeley),
#' lm(cover ~ firesev, data = keeley),
#' lm(firesev ~ age, data = keeley),
#' data = keeley
#' )
#' 
#' update(mod, firesev ~ age + cover)
#' 
#' @export
#' 
update.psem <- function(object, ...) {
  
  l <- list(...)

  for(i in l) {

    if(all(class(i) %in% c("matrix", "data.frame", "SpatialPointsDataFrame", "comparative.data"))) {

      idx <- which(names(object) == "data")

      if(length(idx) == 0) object$data <- i else

        object[[idx]] <- i

      object <- lapply(object, function(j) {

        if(!any(class(j) %in% c("matrix", "data.frame", "SpatialPointsDataFrame", "comparative.data", "formula", "formula.cerror")))

          update(j, data = i) else j

      } )

    } else if(all(class(i) %in% c("character", "formula", "formula.cerror"))) {

      if(length(all.vars.merMod(i)) == 1 | class(i) %in% "formula.cerror") {

        idx <- which(names(object) == "data")

        object <- append(object[-idx], list(i, data = object[[idx]]))

      } else {

        resp <- sapply(object, function(y) if(!any(class(y) %in% c("matrix", "data.frame", "SpatialPointsDataFrame", "comparative.data")))

          all.vars.merMod(y)[1] else "")

        idx <- which(resp == all.vars.merMod(i)[1])

        object[[idx]] <- update(object[[idx]], i)

      }

    } else {

      object[[length(object) + 1]] <- i

    }

  }

  evaluateClasses(object)

  class(object) <- "psem"

  return(object)

}

# New operators for latent, composite variables

# `%~=%`

# `%~+%`
