## Functions from version 1.2 of piecewiseSEM
## For use in vignette()

#' endogenous.reverse
#'
#'  @keywords internal
#'  
endogenous.reverse = function(basis.set, modelList, add.vars = NULL) {

  # Remove NULLs from basis set
  basis.set = basis.set[!sapply(basis.set, is.null)]

  # Identify d-sep tests among endogenous variables and reverse if family != "gaussian"
  names(basis.set) = 1:length(basis.set)

  # Get sorted adjacency matrix
  amat = get.sort.dag(get.formula.list(modelList))

  # Identify intermediate endogenous variables
  idx = colnames(amat)[which(colSums(amat[which(colSums(amat) == 0), , drop = FALSE]) > 0)]

  idx = idx[!idx %in% names(which(colSums(amat[which(!colSums(amat) == 0), , drop = FALSE]) > 0))]

  # Identify variables in the basis set where intermediate endogenous variables are the response
  idx. = sapply(modelList, function(x) all.vars(formula(x))[1] %in% idx)

  # Reverse independence claims
  if(any(idx.)) {

    # Get index of models for which responses are in index
    idx.. = idx[sapply(modelList[idx.], function(x) any(class(x) %in% c("glm", "negbin", "glmmPQL", "glmerMod")))]

    if(length(idx..) > 0) {

      basis.set = append(basis.set,

                         lapply(which(sapply(basis.set, function(i) i[2] %in% idx..)), function(i)

                           c(basis.set[[i]][2], basis.set[[i]][1], basis.set[[i]][-(1:2)])

                         )

      )

    }

  }

  # Ensure that reversal is not attempting to predict exogenous variable
  formulaList = get.formula.list(modelList, add.vars)

  rvars = lapply(formulaList, function(i) {

    if(grepl("cbind\\(.*\\)", paste(formula(i)[2])))

      v = paste0("cbind(", paste(all.vars(i)[1:2], collapse = ","), ")") else

        v = all.vars(i)[1]

      return(gsub(" " , "", v))

  } )

  basis.set = basis.set[sapply(basis.set, function(x) any(x[2] %in% rvars))]

  return(basis.set)

}


#' filter.exogenous
#' 
#' @keywords internal
#' 
filter.exogenous = function(modelList, basis.set, corr.errors = NULL, add.vars = NULL) {

  # Covnert model list into list of vectors
  modelFormulaList = lapply(modelList, function(i) {

    i = formula(i)

    if(grepl("cbind\\(.*\\)", paste(i[2])))

      v = c(paste(formula(i)[2]), all.vars(i)[-(1:2)]) else

        v = all.vars(i)

      return(v)

  } )

  # Get vector of predictor variables
  vars.cleanup = function(x) { x = x[!duplicated(x)]; x = gsub(" ", "" , x); return(x) }

  pred.vars = sapply(modelFormulaList, function(x) x[-1])

  pred.vars = vars.cleanup(pred.vars)

  # Get vector of response variables
  response.vars = sapply(modelFormulaList, function(x) x[1])

  response.vars = vars.cleanup(response.vars)

  # Get vector of variables that appear only as predictors and never as responses
  filter.vars = pred.vars[!pred.vars %in% response.vars]

  # Remove filtered variables when they appear as responses in the basis set
  basis.set = basis.set[!sapply(basis.set, function(i) any(i[2] %in% filter.vars))]

  # Ensure that no entry in the basis set is just the reverse of an existing claim
  basis.set. = lapply(basis.set, function(i) i[order(i)])

  basis.set = basis.set[!duplicated(basis.set.)]

  # Ensure that no entry in the basis set already exists in the model list
  basis.set = lapply(basis.set, function(i)

    if(any(sapply(modelFormulaList, function(j) any(i[1:2] %in% j[1]) & all(i[1:2] %in% j)))) NULL else i

  )

  # Replace interaction : with * when interaction is response in basis set & remove entries attempting to predict interaction
  basis.set = lapply(basis.set, function(i)

    if(is.null(i)) i else { if(i[2] %in% response.vars) gsub(":", "\\*", i) else NULL }

  )

  # Remove NULLs from basis set
  basis.set = basis.set[!sapply(basis.set, is.null)]

  if(length(basis.set) < 1) stop("All endogenous variables are conditionally dependent.\nTest of directed separation not possible!", call. = FALSE)

  return(basis.set)

}

#' get.dag
#' 
#' @keywords internal
#' 
get.dag = function(formulaList) {

  # Insert placeholder for interaction symbol
  formulaList = lapply(formulaList, function(i) formula(gsub("\\:", "_____", Reduce(paste, deparse(i)))))

  # Strip transformations
  formulaList.new = lapply(formulaList, function(i) {

    v = all.vars(i)

    if(grepl("cbind\\(.*\\)", paste(formula(i)[2]))) {

      v = c(paste(v[1], v[2], sep = ","), v[-(1:2)])

    }

    return(v)

  } )

  # Get list of variables and also strip transformations
  vars = unlist(lapply(formulaList, function(i) {

    drop.vars = all.vars(i)[grepl("offset", rownames(attr(terms(i), "factors")))]

    vars. = all.vars(i)[!all.vars(i) %in% drop.vars]

    # Check to see whether response is vector
    if(grepl("cbind\\(.*\\)", paste(formula(i)[2]))) {

      vars. = c(paste(vars.[1], vars.[2], sep = ","), vars.[-(1:2)])

    }

    return(vars.)

  } ) )

  vars = unname(vars[!duplicated(vars)])

  # Create adjacency matrix
  amat = do.call(cbind, lapply(vars, function(i) {

    # Identify formula where variable is response
    form.no = which(sapply(formulaList.new, function(j) j[1] == i))

    if(length(form.no) == 0) rep(0, length(vars)) else {

      # Isolate variable from formula list
      form = formulaList.new[[which(sapply(formulaList.new, function(j)j[1] == i))]]

      vars %in% form + 0

    }

  } ) )

  # Ensure diagonal is zero
  diag(amat) = 0

  # Name rows and columns
  dimnames(amat) = list(vars, vars)

  # Determine if graph is acylic
  if(all(colSums(amat) > 0)) stop("Model is non-recursive (cyclical)! Remove directed cycles and re-run.", call. = FALSE)

  return(amat)

}

#' get.basis.set
#' 
#' @keywords internal
#' 
get.basis.set = function(amat) {

  ret = lapply(1:ncol(amat), function(i) {

    lapply(i:nrow(amat), function(j) {

      if(amat[i, j] != 0 | i == j) NULL else {

        # Get variables for independence test
        dsep = unlist(dimnames(amat[i, j, drop = FALSE]))

        # Get vector of conditional variables
        cond.var = c(
          rownames(amat)[which(amat[, dsep[1], drop = FALSE] == 1)],
          rownames(amat)[which(amat[, dsep[2], drop = FALSE] == 1)]
        )

        # Remove conditional variables already in the independence claim
        cond.var = cond.var[!cond.var %in% dsep]

        # Return full independence claim
        c(dsep, cond.var)

      }

    } )

  } )

  ret = unlist(ret, recursive = FALSE)

  ret = lapply(ret, function(j) j[!duplicated(j)])

  ret = ret[!sapply(ret, is.null)]

  if(length(ret) == 0)

    stop("All endogenous variables are conditionally dependent.\nTests of directed separation not possible!", call. = FALSE)

  # Add binding for binomial variables
  for(i in 1:length(ret)) {

    if(any(grepl(",", ret[[i]]))) {

      idx = which(grepl(",", ret[[i]]))

      for(j in idx) {

        ret[[i]][j] = paste0("cbind(", ret[[i]][j], ")")

      }

    }

  }

  return(ret)

}

#' get.formula.list
#' 
#' @keywords internal
#' 
get.formula.list = function(modelList, add.vars = NULL) {

  # Get list of formula from model list
  formulaList = lapply(modelList, function(i)

    if(all(class(i) %in% c("lm", "rq", "glm", "negbin", "lme", "glmmPQL", "gls", "pgls", "glmmadmb"))) formula(i) else

      if(all(class(i) %in% c("lmerMod", "merModLmerTest", "glmerMod", "glmmTMB"))) nobars(formula(i))

  )

  if(any(unlist(lapply(formulaList, is.null)))) stop("At least one model class not yet supported", .call = FALSE)

  # Check to see if any variables in the formula list appear in add.vars
  if(any(unlist(lapply(formulaList, all.vars)) %in% add.vars)) stop("Variables in the model list appear in add.vars!")

  # If additional variables are present, add them to the basis set
  if(!is.null(add.vars)) {

    formulaList = append(formulaList, unname(sapply(add.vars, function(x) as.formula(paste(x, x, sep = "~")))))

  }

  # Expand interactions in formula list
  formulaList = lapply(formulaList, function(i)

    if(grepl("\\*|:", paste(format(formula(i)), collapse = ""))) {

      lhs = paste(rownames(attr(terms(i), "factors"))[1])

      rhs = attr(terms(i), "term.labels")

      # Sort interactions so always alphabetical
      for(j in which(grepl(":", rhs))) {

        # Split interactions and sort alphabetically
        int = unlist(lapply(strsplit(rhs[j], ":"), sort))

        # Recombine
        int.rec = paste(int, collapse = ":")

        # Re-insert into formula
        rhs[j] = int.rec

      }

      # Collapse into formula
      rhs = paste(rhs, collapse = " + ")

      # And return full formula
      formula(paste(lhs, " ~ ", rhs))

    }

    else i

  )

  if(any(sapply(formulaList, function(x) grepl("poly", x))))

    stop("Polynomials computed within the regression are not yet allowed.\nCompute externally and supply each component as a separate variable!")

  return(formulaList)

}

#' get.model.control
#' 
#' @keywords internal
#' 
get.model.control = function(model, model.control) {

  model.class = if(inherits(model, "merModLmerTest")) "lmerMod" else class(model)

  # Match model control list with appropriate model class for basis model
  if(is.null(model.control)) {

    if(inherits(model, "glm")) glm.control()

    else if(inherits(model, "gls")) glsControl()

    else if(any(class(model) %in% c("lme", "glmmPQL", "glmmadmb"))) lmeControl()

    else if(any(class(model) %in% c("lmerMod", "merModLmerTest"))) lmerControl()

    else if(inherits(model, "glmerMod")) glmerControl()

  } else {

    ## assume model.control = a *list* of control lists
    control.classes = lapply(model.control, function(i)

      gsub("(.*)Control", "\\1", class(i))[1] )

    if(any(control.classes %in% (M.cl <- gsub("(.*)Mod", "\\1", model.class))))

      model.control[[control.classes %in% M.cl]] else

        if(inherits(model, "glm"))

          model.control[[sapply(model.control, length) >= 3]] else

            if(inherits(model, "gls"))

              model.control[[sapply(model.control, length) >= 13]] else

                if(any(class(model) %in% c("lme", "glmmPQL", "glmmadmb")))

                  model.control[[sapply(model.control, length) >= 15]]

  }

}

#' get.random.formula
#' 
#' @keywords internal
#' 
get.random.formula = function(model, rhs, modelList, dropterms = NULL) {

  if(class(rhs) == "formula") rhs = Reduce(paste, deparse(rhs))

  # Get random formula from model
  random.formula = if(any(class(model) %in% c("lme", "glmmPQL", "glmmadmb")))

    deparse(model$call$random) else

      if(any(class(model) %in% c("lmerMod", "merModLmerTest", "glmerMod", "glmmTMB")))

        paste("(", findbars(formula(model)), ")", collapse = " + ")

  # Get random structure(s)
  random.structure = if(any(class(model) %in% c("lme", "glmmPQL", "glmmadmb"))) {

    # If crossed random effects, extract random effects and store in a list
    if(grepl("list\\(", as.character(random.formula)))

      as.list(
        gsub(
          "(.*)=.*",
          "\\1",
          strsplit(gsub("list\\((.*)\\)", "\\1", random.formula), ",")[[1]]
        )
      )

    else gsub(" ", "", gsub(".*\\|(.*)?", "\\1", random.formula))

  } else

    if(any(class(model) %in% c("lmerMod", "merModLmerTest", "glmerMod", "glmmTMB"))) {

      ran.ef.splt = strsplit(random.formula, "\\+.")[[1]]

      sapply(ran.ef.splt[sapply(ran.ef.splt, function(x) grepl("\\|", x))],

             function(x)

               gsub(" ", "", gsub(".*\\|(.*)\\)", "\\1", x))

      )

    }

  random.structure = unname(random.structure[!duplicated(random.structure)])

  # Get random slopes in the model list, otherwise return vector of terms to drop
  random.slopes =

    if(any(class(model) %in% c("lme", "glmmPQL", "glmerMod", "merModLmerTest", "glmmadmb", "glmmTMB")))

      if(is.null(dropterms)) {

        unlist(lapply(1:length(modelList), function(i) {

          if(any(class(modelList[[i]]) %in% c("glmmPQL"))) {

            ran.ef = ifelse(any(class(modelList[[i]]$coefficients$random) != "list"),

                            list(modelList[[i]]$coefficients$random),

                            modelList[[i]]$coefficients$random)

            as.vector(sapply(ran.ef, function(j) colnames(j)))

          }

          else if(any(class(modelList[[i]]) %in% c("lme", "glmerMod", "merModLmerTest", "glmmadmb", "glmmTMB"))) {

            ran.ef = ifelse(any(class(ranef(modelList[[i]])) != "list"),

                            list(ranef(modelList[[i]])),

                            ranef(modelList[[i]]))

            as.vector(sapply(ran.ef, function(j) colnames(j)))

          }

        } ) )

      } else dropterms

  random.slopes = unname(random.slopes[!duplicated(random.slopes) & random.slopes != "(Intercept)"])

  # Define new random slopes
  new.random.slopes = random.slopes[which(random.slopes %in% unlist(strsplit(rhs, ".\\+.")))]

  if(length(new.random.slopes) == 0) new.random.slopes = 1 else new.random.slopes = paste0(new.random.slopes, collapse = " + ")

  # Replace random slopes if any variables in model formula appear in random slopes
  if(length(random.slopes) != 0) {

    if(any(class(model) %in% c("lme", "glmmPQL", "glmmadmb")))

      if(is.list(random.structure)) {

        eval(parse(text = gsub("*\\~(.*)", paste0("~ ", new.random.slopes, "))"), random.formula)))

      } else {

        formula(
          paste("~ ",
                new.random.slopes,
                " | ",
                random.structure)
        )

      } else if(any(class(model) %in% c("glmerMod", "merModLmerTest", "glmmTMB"))) {

        paste(
          sapply(random.structure, function(x)
            paste("(", new.random.slopes, " | ", x, ")") ),
          collapse = " + ")

      }

  } else if(length(random.slopes) == 0) {

    if(is.list(random.structure)) {

      eval(parse(text = gsub("*\\~(.*)", paste0("~ ", new.random.slopes, "))"), random.formula)))

    } else if(any(class(model) %in% c("lme", "glmmPQL"))) formula(random.formula)

    else random.formula

  }

}

#' get.scaled.data
#' 
#' @keywords internal
#' 
get.scaled.data = function(modelList, data, standardize) {

  if(is.null(data)) stop("Must supply data if calculating standardized coefficients!")

  if(any(sapply(modelList, class) == "pgls") | class(data) == "comparative.data") {

    # Extract data.frame
    newdata = data$data

  } else {

    newdata = data

  }

  # Identify variables that are transformed in the model formula
  transform.vars = unlist(lapply(modelList, function(i) {

    # Break apart formula
    if(any(class(i) == "pgls")) vars = i$varNames else

      vars = rownames(attr(terms(i), "factors"))

    # Identify transformed vars
    vars[grepl("log\\(|log10\\(|sqrt\\(|I\\(", vars)]

  } ) )

  # Remove duplicates
  transform.vars = transform.vars[!duplicated(transform.vars)]

  # Strip transformations
  transform.vars2 = sapply(transform.vars, function(i)

    gsub(" ", "",
         gsub(".*([[:alpha:]]).*", "\\1",
              gsub(".*\\((.*)\\).*", "\\1", i)
         )
    )

  )

  # For each variables in transform.vars, perform transformation and store as original variable
  if(length(transform.vars) > 0)

    for(i in 1:length(transform.vars2)) {

      # Perform transformation
      newdata[, transform.vars2[i]] =

        sapply(newdata[, transform.vars2[i]], function(x)

          eval(parse(text = gsub(transform.vars2[i], x, transform.vars[i])))

        )

    }

  # Get variables to scale, ignoring variables that are modeled to non-normal distributions
  vars = unlist(lapply(modelList, function(x) {

    if(grepl("cbind", deparse(formula(x))))

      all.vars(formula(x))[-c(1:2)] else

        all.vars(formula(x))

  } ) )

  vars = vars[!duplicated(vars)]

  non.normal.vars = unlist(lapply(modelList, function(i) {

    family = if(any(class(i) %in% c("glmmadmb"))) i$family else

      if(any(class(i) %in% c("glm", "glmerMod", "negbin"))) family(i) else

        NULL

    if(!is.null(family)) all.vars(formula(i))[1]

  } ) )

  vars.to.scale = vars[!vars %in% non.normal.vars]

  if(!is.null(non.normal.vars))

    warning("One or more responses not modeled to a normal distribution: keeping response(s) on original scale!")

  # Remove variables that are factors
  vars.to.scale = vars.to.scale[!vars.to.scale %in% colnames(data)[sapply(data, function(x) any(is.factor(x) | is.character(x)))] ]

  # Convert variables that are factors to numeric
  # data[, sapply(data, function(x) any(is.factor(x) | is.character(x)))] =
  #
  #   apply(data[, sapply(data, function(x) any(is.factor(x) | is.character(x))), drop = FALSE], 2, function(x) as.numeric(as.factor(x)))

  # Remove variables that appear as random effects
  rand.mods = which(sapply(modelList, class) %in% c("lme", "lmerMod", "merModLmerTest", "glmerMod"))

  rand.effs = c()

  for(i in rand.mods) rand.effs = c(rand.effs, names(ranef(modelList[[i]])))

  # Unnest nested variables
  if(length(rand.effs) > 0) {

    rand.effs = unlist(strsplit(rand.effs, ":"))

    vars.to.scale = vars.to.scale[!vars.to.scale %in% rand.effs]

  }

  # Remove duplicated variables
  vars.to.scale = vars.to.scale[!duplicated(vars.to.scale)]

  # Remove transformed variables already scaled
  vars.to.scale = vars.to.scale[!vars.to.scale %in% gsub(" ", "", transform.vars)]

  # Run check to see if variables appear as columns
  if(!all(vars.to.scale %in% colnames(newdata))) stop("Some predictors do not appear in the dataset!")

  # Scale those variables by mean and SD, or by range
  newdata[, vars.to.scale] = apply(newdata[, vars.to.scale, drop = FALSE], 2, function(x) {

    if(standardize == "scale") scale(x) else

      if(standardize == "range") (x-min(x, na.rm = T)) / diff(range(x, na.rm = T)) else

        x

  } )

  if(class(data) == "comparative.data") {

    data$data = newdata

    newdata = data

  }

  return(newdata)

}

#' get.scaled.model
#' 
#' @keywords internal
#' 
get.scaled.model = function(model, newdata, modelList) {

  if(any(class(model) %in% c("lmerMod", "merModLmerTest", "glmerMod", "glmmTMB"))) {

    # Get random effects
    rand.effs = gsub(" ", "", sapply(findbars(formula(model)), function(x) gsub(".*\\|(.*)", "\\1", deparse(x))))

    # Unnest nested variables
    rand.effs = unlist(strsplit(rand.effs, ":"))

    # Get fixed effects
    fixed.effs = all.vars(formula(model))[!all.vars(formula(model)) %in% rand.effs]

    # Get fixed formula stripped of transformations
    if(grepl("cbind", deparse(formula(model))))

      fixed.form = paste0("cbind(", fixed.effs[1], ", ", fixed.effs[2], ") ~ ", paste0(fixed.effs[-(1:2)], collapse = " + ")) else

        fixed.form = paste0(fixed.effs[1], " ~ ", paste0(fixed.effs[-1], collapse = " + "))

    # Bind back in random structure
    random.form = get.random.formula(model, rhs = paste0(fixed.effs[-1], collapse = " + "), modelList)

    # Get updated formula
    new.form = paste0(fixed.form, " + ", random.form)

  } else {

    new.form = paste(all.vars(formula(model))[1], paste(formula(model)[-2], collapse = ""))

  }

  # Update model
  if(any(class(model) %in% c("lme", "glmmPQL")))

    model = update(model, fixed = as.formula(new.form), data = newdata) else

      model = update(model, as.formula(new.form), data = newdata)

    # Return model
    return(model)

}

#' get.sort.dag
#' 
#' @keywords internal
get.sort.dag = function(formulaList) {

  # Get adjaceny matrix
  amat = get.dag(formulaList)

  # Get predictors where colSums == 0
  idx = unname(which(colSums(amat) == 0))

  # Of remaining variables, look for only those with links to the predictors in col.zero
  idx. =
    which(!colnames(amat) %in% colnames(amat)[idx])[
      colSums(amat[!colnames(amat) %in% colnames(amat)[idx], !colnames(amat) %in% colnames(amat)[idx], drop = FALSE]) == 0 &
        colSums(amat[colnames(amat) %in% colnames(amat)[idx], !colnames(amat) %in% colnames(amat)[idx], drop = FALSE]) != 0]

  # Sort by increasing position in amat
  pos = apply(amat[, idx., drop = FALSE], 2, function(x) max(which(x == 1)))

  idx. = c(idx, which(colnames(amat) %in% names(pos[order(pos)])))

  # Sort remaining variables by increasing position in amat
  pos. = apply(amat[, !colnames(amat) %in% colnames(amat)[idx.], drop = FALSE], 2, function(x) max(which(x == 1)))

  idx.. = c(idx., which(colnames(amat) %in% names(pos.[order(pos.)])))

  # Return adjacency matrix, sorted
  amat[idx.., idx..]

}

#' partial.resid
#' 
#' @keywords internal
#' 
partial.resid <- function(...) {

  warning("`partial.resid` has been replaced by `partialResid`", call. = FALSE)

}

#' sem.aic
#' 
#' @keywords internal
#' 
sem.aic = function(

  modelList, data, corr.errors = NULL, add.vars = NULL, grouping.vars = NULL, grouping.fun = mean,
  adjust.p = FALSE, basis.set = NULL, pvalues.df = NULL, model.control = NULL, .progressBar = TRUE

) {

  if(is.null(basis.set)) basis.set = suppressWarnings(sem.basis.set(modelList, corr.errors, add.vars))

  if(is.null(pvalues.df)) pvalues.df = suppressMessages(suppressWarnings(sem.missing.paths(

    modelList, data, conditional = FALSE, corr.errors, add.vars, grouping.vars,
    grouping.fun, adjust.p, basis.set, model.control, .progressBar

  ) ) )

  fisher.c = sem.fisher.c(

    modelList, data, corr.errors, add.vars, grouping.vars, grouping.fun,
    adjust.p, basis.set, pvalues.df, model.control, .progressBar

  )

  # Calculate likelihood degrees of freedom
  K = do.call(sum, lapply(modelList, function(i) attr(logLik(i), "df")))

  # Calculate AIC
  AIC = unname(fisher.c[1] + 2 * K)

  # Calculate AICc
  n.obs = min(sapply(modelList, function(x) {

    if(class(x) == "rq")

      length(na.omit(residuals(x))) else

        nobs(x)

  } ) )

  AICc = unname(fisher.c[1] + 2 * K * (n.obs/(n.obs - K - 1)))

  # Return output in a data.frame
  data.frame(
    AIC = round(AIC, 3),
    AICc = round(AICc, 3),
    K = round(K, 1),
    n = round(n.obs, 1) )

}

#' sem.basis.set
#' 
#' @keywords internal
#' 
sem.basis.set = function(modelList, corr.errors = NULL, add.vars = NULL) {

  # Get list of formula from model list
  formulaList = get.formula.list(modelList, add.vars)

  # Stop if any response in the model list is the same as any other response
  if(any(duplicated(sapply(formulaList, function(x) all.vars(x)[1]))))

    stop("Duplicate responses detected in the model list.\n
         Collapse multiple single regressions into a single multiple regression so that each response appears only once!")

  # Get adjacency matrix and sort by parent to child nodes
  amat = get.sort.dag(formulaList)

  # If intercept only model, add response variable to adjacency matrix
  if(any(unlist(lapply(modelList, function(i) deparse(formula(i)[2]) %in% c("~1", "~ 1"))))) {

    # Isolate intercept only model(s)
    responses = sapply(modelList[which(sapply(modelList, function(i) grepl("~ 1|~1", deparse(formula(i)))))],

                       function(j) strsplit(paste(formula(j)), "~")[[2]]

    )

    amat = cbind(

      rbind(amat, matrix(rep(0, dim(amat)[1]), nrow = 1, dimnames = list(responses))),

      matrix(rep(0, dim(amat)[1] + 1), ncol = 1, dimnames = list(NULL, responses))

    )

  }

  # Generate basis set
  if(all(amat == 0) & all(dim(amat) == 1)) basis.set = NULL else basis.set = get.basis.set(amat)

  # Replace placeholder for interaction symbol with :
  basis.set = lapply(basis.set, function(i) gsub(paste("_____", collapse = ""), "\\:", i))

  # If correlated errors are present, remove them from the basis set
  if(!is.null(corr.errors)) {

    basis.set =  lapply(1:length(basis.set), function(i) {

      inset = unlist(lapply(corr.errors, function(j) {

        corr.vars = gsub(" ", "", unlist(strsplit(j,"~~")))

        all(

          unlist(

            lapply(1:2, function(k)

              grepl(paste(corr.vars, collapse = "|"), basis.set[[i]][k])

            )
          )
        )

      } ) )

      if(any(inset == TRUE)) NULL else basis.set[[i]]

    } )
  }

  # Replace any d-sep where interactions are regressed against the main effect with NULL
  basis.set = lapply(basis.set, function(i) {

    if(is.null(i)) NULL else {

      if(grepl("\\*|\\:", i[1])) {

        int = strsplit(i[1], "\\*|\\:")[[1]]

        if(any(int %in% i[2])) NULL else i

      }

      else i

    }

  } )

  # Identify responses for which offset is present
  rpl = do.call(rbind, lapply(formulaList, function(i) {

    lhs = paste(rownames(attr(terms(i), "factors"))[1])

    rhs = rownames(attr(terms(i), "factors"))[-1]

    if(any(grepl("offset", rhs)))

      data.frame(response = lhs, offset = rhs[grepl("offset", rhs)]) else

        NULL

  } ) )

  # Add offset to basis set
  if(!is.null(rpl))

    basis.set = lapply(basis.set, function(i) {

      if(any(i[2] == rpl$response)) {

        c(i, as.character(rpl[rpl$response == i[2], "offset"]))

      } else i

    } )


  ### TEMPORARY FIX ###

  # Reverse intermediate endogenous variables fitted to non-normal distributions
  basis.set = endogenous.reverse(basis.set, modelList, add.vars)

  ### TEMPORARY FIX ###


  # Filter exogenous predictors from the basis set
  basis.set = filter.exogenous(modelList, basis.set, corr.errors, add.vars)

  # Re-apply transformations
  basis.set = lapply(basis.set, function(i) {

    # Get list of transformed predictors
    t.pvars = lapply(formulaList, function(x) colnames(attr(terms(x), "factors")))

    # Get list of untransformed predictors
    pvars = lapply(formulaList, function(i) {

      if(grepl("cbind\\(.*\\)", paste(formula(i)[2])))

        v = all.vars(i)[-(1:2)] else

          v = all.vars(i)[-1]

        return(v)

    } )

    # Get list of transformed responses
    t.rvars = lapply(formulaList, function(x) gsub(" " , "", rownames(attr(terms(x), "factors"))[1]))

    # Get list of untransformed responses
    rvars = lapply(formulaList, function(i) {

      if(grepl("cbind\\(.*\\)", paste(formula(i)[2])))

        v = paste0("cbind(", paste(all.vars(i)[1:2], collapse = ","), ")") else

          v = all.vars(i)[1]

        return(gsub(" " , "", v))

    } )

    # Re-transform predictors
    for(j in (1:length(i))[-2]) {

      # Get variable index for lookup
      idx = cbind(

        which(sapply(pvars, function(k) any(k %in% i[j]))),

        unlist(sapply(pvars, function(l) which(l == i[j])))

      )

      if(sum(idx) > 0) {

        # Conduct lookup
        t.pvar = sapply(1:nrow(idx), function(m) t.pvars[[idx[m, 1]]][idx[m, 2]] )

        t.pvar = t.pvar[!duplicated(t.pvar)]

        # Replace predictors
        if(length(t.pvar) > 0) i[j] = t.pvar[which.max(sapply(t.pvar, function(p) nchar(p)))]

      }

    }

    # Get variable index for lookup
    idx = cbind(

      which(sapply(rvars, function(k) any(k %in% i[2]))),

      unlist(sapply(rvars, function(l) which(l == i[2])))

    )

    if(sum(idx) != 0) {

      # Conduct lookup
      t.rvar = sapply(1:nrow(idx), function(m) t.rvars[[idx[m, 1]]][idx[m, 2]] )

      t.rvar = t.rvar[!duplicated(t.rvar)]

      # Replace predictors
      if(length(t.rvar) > 0) i[2] = t.rvar[which.max(sapply(t.rvar, function(p) nchar(p)))]

    }

    return(i)

  } )

  # Remove NULLs from basis set
  basis.set = basis.set[!sapply(basis.set, is.null)]

  # Re-assign names from dropped entries
  names(basis.set) = as.numeric(as.factor(as.numeric(names(basis.set))))

  return(basis.set)

}

#' sem.fisher.c
#' 
#' @keywords internal
#' 
sem.fisher.c = function(

  modelList, data, corr.errors = NULL, add.vars = NULL, grouping.vars = NULL, grouping.fun = mean,
  adjust.p = FALSE, basis.set = NULL, pvalues.df = NULL, model.control = NULL, .progressBar = TRUE

) {

  if(is.null(basis.set)) basis.set = suppressWarnings(sem.basis.set(modelList, corr.errors, add.vars))

  if(is.null(pvalues.df)) pvalues.df = suppressMessages(suppressWarnings(sem.missing.paths(

    modelList, data, conditional = FALSE, corr.errors, add.vars, grouping.vars,
    grouping.fun, adjust.p, basis.set, model.control, .progressBar

  ) ) )

  # Convert any p-values to a very small number as log(0) == -Inf
  if(length(basis.set) > 0 & any(pvalues.df$p.value == 0)) pvalues.df[pvalues.df$p.value == 0, "p.value"] = 2e-16

  # Calculate Fisher's C statistic
  fisher.C = if(length(basis.set) > 0) -2 * sum(log(pvalues.df$p.value)) else 0
  # Calculate associated p-value from Chi-squared distribution
  p.value = 1 - pchisq(fisher.C, 2 * length(basis.set))

  # Return output in a data.frame
  data.frame(fisher.c = round(fisher.C, 2), df = round(2 * length(basis.set), 1), p.value = round(p.value, 3))

}

#' sem.coefs
#' 
#' @keywords internal
#' 
sem.coefs = function(modelList, data = NULL, standardize = "none", corr.errors = NULL, intercept = FALSE) {

  warning("`sem.coefs` has been replaced. Use `psem` instead of `list`, and then call `summary` or `coefs` on that object", call. = FALSE)

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

    } else if(any(class(i) %in% c("rq"))) {

      tab = summary(i, se = "boot")$coefficients

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

    } else if(any(class(i) %in% c("glmmTMB"))) {

      tab = summary(i)$coefficients$cond

      data.frame(response = Reduce(paste, deparse(formula(i)[[2]])),
                 predictor = rownames(tab)[irow],
                 estimate = tab[irow, 1],
                 std.error = tab[irow, 2],
                 p.value = tab[irow, 4]
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

#' sem.fit
#' 
#' @keywords internal
#' 
sem.fit = function(

  modelList, data, conditional = FALSE, corr.errors = NULL, add.vars = NULL, grouping.vars = NULL,
  grouping.fun = mean, adjust.p = FALSE, basis.set = NULL, pvalues.df = NULL, model.control = NULL,
  .progressBar = TRUE

) {

  warning("`sem.fit` has been replaced. Use `psem` instead of `list`, and then call `summary` on that object", call. = FALSE)

  if(is.null(data)) stop("Must supply dataset")

  n.obs = sapply(modelList, function(x) {

    if(class(x) == "rq")

      length(na.omit(residuals(x))) else

        nobs(x)

  } )

  if(!all(n.obs)) warning("All models do not have the same number of observations")

  # Get basis set
  if(is.null(basis.set)) basis.set = suppressMessages(suppressWarnings(

    sem.basis.set(modelList, corr.errors, add.vars)

  ) )

  # Conduct d-sep tests
  if(is.null(pvalues.df)) pvalues.df = sem.missing.paths(

    modelList, data, conditional, corr.errors, add.vars, grouping.vars,
    grouping.fun, adjust.p, basis.set, model.control, .progressBar

  )

  # Derive Fisher's C statistic and compare to Chi-squared distribution
  fisher.c = sem.fisher.c(

    modelList, data, corr.errors, add.vars, grouping.vars,
    grouping.fun, adjust.p, basis.set, pvalues.df, model.control, .progressBar

  )

  # Use Fisher's C to derive AIC values
  AIC.c = sem.aic(

    modelList, data, corr.errors, add.vars, grouping.vars,
    grouping.fun, adjust.p, basis.set, pvalues.df, model.control, .progressBar

  )

  # Round values in output table
  pvalues.df[, c(2:3, 5:6)] = apply(pvalues.df[, c(2:3, 5:6)], 2, function(x) round(x, 4) )

  # Return d-sep tests, Fisher's C, and AIC values in a list
  l = list(pvalues.df, fisher.c, AIC.c)

  names(l) = c("missing.paths", "Fisher.C", "AIC")

  return(l)

}

#' sem.missing.paths
#' 
#' @keywords internal
#' 
sem.missing.paths = function(

  modelList, data, conditional = FALSE, corr.errors = NULL, add.vars = NULL, grouping.vars = NULL,
  grouping.fun = mean, adjust.p = FALSE, basis.set = NULL, model.control = NULL, .progressBar = TRUE

) {

  # Get basis set
  if(is.null(basis.set)) basis.set = suppressWarnings(sem.basis.set(modelList, corr.errors, add.vars))

  # Add progress bar
  if(.progressBar == T & length(basis.set) > 0) pb = txtProgressBar(min = 0, max = length(basis.set), style = 3) else pb = NULL

  # Perform d-sep tests
  if(length(basis.set) > 0) pvalues.df = do.call(rbind, lapply(1:length(basis.set), function(i) {

    # Get basis model from which to build the d-sep test
    basis.mod = modelList[[which(sapply(modelList, function(j) {

      if(any(class(j) == "pgls")) j = j$formula

      gsub(" ", "", rownames(attr(terms(j), "factors"))[1]) == basis.set[[i]][2]

    } ) ) ]]

    # Get fixed formula
    rhs = if(length(basis.set[[i]]) <= 2) paste(basis.set[[i]][1]) else

      paste(basis.set[[i]][c(3:length(basis.set[[i]]), 1)], collapse = " + ")

    # Get random formula
    random.formula = get.random.formula(basis.mod, rhs, modelList)

    # Aggregate at the level of the grouping variable
    if(!is.null(grouping.vars)) {

      # Test to see if response is identical for levels of grouping factor
      response = as.character(formula(basis.mod)[2])

      response.test = by(data, lapply(grouping.vars, function(j) data[ ,j]), function(x) length(unique(x[, response])) == 1)

      response.test = na.omit(vapply(response.test, unlist, unlist(response.test[[1]])))

      # If so, aggregate by grouping.vars
      if(all(response.test == TRUE)) {

        # Get named list
        groups = lapply(grouping.vars, function(j) data[ ,j])

        names(groups) = grouping.vars

        # Decide aggregation level
        groups = groups[which(sapply(groups, function(x) length(unique(x))) == length(response.test)):length(groups)]

        # Aggregate and replace data
        data = suppressWarnings(

          aggregate(data, by = groups, grouping.fun, na.rm = T)

        )

        # Remove duplicated colnames (keep first instance, from above)
        data = data[, !duplicated(colnames(data))]

      }

    }

    # Get model controls
    control = get.model.control(basis.mod, model.control)

    # Update basis model with new formula and random structure based on d-sep
    basis.mod.new = suppressMessages(suppressWarnings(

      if(is.null(random.formula) & class(basis.mod) == "rq")

        update(basis.mod, formula(paste(basis.set[[i]][2], " ~ ", rhs)), data = data) else

          if(is.null(random.formula) | class(basis.mod) == "glmmadmb")

            update(basis.mod, formula(paste(basis.set[[i]][2], " ~ ", rhs)), control = control, data = data) else

              if(any(class(basis.mod) %in% c("lme", "glmmPQL")))

                update(basis.mod, fixed = formula(paste(basis.set[[i]][2], " ~ ", rhs)), random = random.formula, control = control, data = data) else

                  if(any(class(basis.mod) %in% "glmmTMB")) update(basis.mod, formula = formula(paste(basis.set[[i]][2], " ~ ", rhs, " + ", random.formula)), data = data) else

                    update(basis.mod, formula = formula(paste(basis.set[[i]][2], " ~ ", rhs, " + ", random.formula)), control = control, data = data)

    ) )

    # Get row number from coefficient table for d-sep variable
    # if(any(!class(basis.mod.new) %in% c("pgls"))) {
    #
    #   # Get row number of d-sep claim
    #   row.num = which(basis.set[[i]][1] == rownames(attr(terms(basis.mod.new), "factors"))[-1]) + 1
    #
    #   # Get row number if interaction variables are switched
    #   if(length(row.num) == 0 & grepl("\\:|\\*", basis.set[[i]][1])) {
    #
    #     # If interaction is reported as asterisk, convert to semicolon
    #     int = gsub(" \\* ", "\\:", basis.set[[i]][1])
    #
    #     # Get all combinations of interactions
    #     all.ints = sapply(strsplit(int, ":"), function(x) {
    #
    #       datf = expand.grid(rep(list(x), length(x)), stringsAsFactors = FALSE)
    #
    #       datf = datf[apply(datf, 1, function(x) !any(duplicated(x))), ]
    #
    #       apply(datf, 1, function(x) paste(x, collapse = ":"))
    #
    #     } )
    #
    #     row.num = which(attr(terms(basis.mod.new), "term.labels") %in% all.ints) + 1
    #
    #     }
    #
    #   } else {
    #
    #     row.num = which(basis.set[[i]][1] == basis.mod.new$varNames)
    #
    #   }

    # Return new coefficient table
    ret = if(any(class(basis.mod.new) %in% c("lmerMod", "merModLmerTest"))) {

      coef.table = suppressMessages(summary(basis.mod.new)$coefficients)

      # Get P-values baesd on Kenward-Rogers approximation of denominator degrees of freedom
      basis.mod.drop = update(basis.mod.new, as.formula(paste("~ . -", basis.set[[i]][1])))

      kr.p = KRmodcomp(basis.mod.new, basis.mod.drop)

      # Combine with coefficients from regular ouput
      data.frame(
        t(coef.table[nrow(coef.table), 1:2]),
        kr.p$test$ddf[1],
        coef.table[nrow(coef.table), 3],
        kr.p$test$p.value[1],
        row.names = NULL
      )

    } else if(any(class(basis.mod.new) %in% c("lm", "glm", "negbin", "pgls", "glmerMod", "glmmadmb"))) {

      coef.table = summary(basis.mod.new)$coefficients

    } else if(any(class(basis.mod.new) %in% c("rq"))) {

      coef.table = summary(basis.mod.new, se = "boot")$coefficients

    } else if(any(class(basis.mod.new) == "glmmTMB")) {

      coef.table = summary(basis.mod.new)$coefficients$cond

    } else {

      coef.table = summary(basis.mod.new)$tTable

    }

    if(!any(class(basis.mod.new) %in% c("lmerMod", "merModLmerTest")))

      ret = as.data.frame(t(unname(coef.table[nrow(coef.table), ])))

    # Add df if summary table does not return
    if(length(ret) != 5 & any(class(basis.mod.new) %in% c("lm", "glm", "negbin", "pgls")))

      ret = cbind(ret[1:2], summary(basis.mod.new)$df[2], ret[3:4]) else

        if(length(ret) != 5 & any(class(basis.mod.new) %in% c("rq")))

          ret = cbind(ret[1:2], summary(basis.mod.new, se = "boot")$rdf, ret[3:4]) else

            if(length(ret) != 5 & any(class(basis.mod.new) %in% c("glmmadmb")))

              ret = cbind(ret[1:2], summary(basis.mod.new)$n, ret[3:4]) else

                if(length(ret) != 5)

                  ret = cbind(ret[1:2], NA, ret[3:4])

    # Rename columns
    names(ret) = c("estimate", "std.error", "df", "crit.value", "p.value")

    # Adjust p-value based on Shipley 2013
    if(adjust.p == TRUE) {

      if(any(class(basis.mod.new) %in% c("lme", "glmmPQL"))) {

        t.value = coef.table[nrow(coef.table), 4]

        ret[5] = 2*(1 - pt(abs(t.value), nobs(basis.mod.new) - sum(apply(basis.mod.new$groups, 2, function(x) length(unique(x))))))

      } else if(any(class(basis.mod.new) %in% c("lmerMod", "glmerMod"))) {

        z.value = coef.table[nrow(coef.table), 3]

        ret[5] = 2*(1 - pt(abs(z.value), nobs(basis.mod.new) - sum(summary(basis.mod.new)$ngrps)))

      } else if(any(class(basis.mod.new) %in% c("glmmadmb"))) {

        z.value = coef.table[nrow(coef.table), 3]

        ret[5] = 2*(1 - pt(abs(z.value), nobs(basis.mod.new) - sum(summary(basis.mod.new)$npar)))

      } else if(any(class(basis.mod.new) %in% c("glmmTMB"))) {

        z.value = coef.table[nrow(coef.table), 3]

        ret[5] = 2*(1 - pt(abs(z.value), nobs(basis.mod.new) - sum(summary(basis.mod.new)$ngrps$cond)))

      }

    }

    if(.progressBar == TRUE) setTxtProgressBar(pb, i)

    # Modify rhs if number of characters exceeds 20
    rhs.new =

      if(length(basis.set[[i]]) < 3) rhs else {

        if(conditional == FALSE)

          paste0(basis.set[[i]][1], " + ...") else

            paste(basis.set[[i]][c(1, 3:length(basis.set[[i]]))], collapse = " + ")

      }

    # Bind in d-sep metadata
    data.frame(missing.path = paste(basis.set[[i]][2], " ~ ", rhs.new, sep = ""), ret)

  } ) ) else

    pvalues.df = data.frame(missing.path = NA, estimate = NA, std.error = NA, DF = NA, crit.value = NA, p.value = NA)

  # Identify duplicate tests from intermediate endogenous variables
  dup = names(basis.set)

  # Return lowest P-value
  if(any(duplicated(dup))) {

    pvalues.df = do.call(rbind, lapply(unique(dup), function(x) {

      if(length(dup[dup == as.numeric(x)]) > 1)

        warning("Some d-sep tests are non-symmetrical. The most conservative P-value has been returned. Stay tuned for future developments...")

      pvalues.df[as.numeric(x), ][which.min(pvalues.df[as.numeric(x), "p.value"]), ]

    }

    ) )

  }

  # Set degrees of freedom as numeric
  pvalues.df$df = round(as.numeric(pvalues.df$df), 1)

  if(!is.null(pb)) close(pb)

  if(any(grepl("...", pvalues.df$missing.path)) & conditional != TRUE)

    message("Conditional variables have been omitted from output table for clarity (or use argument conditional = T)")

  # rm(dup)

  # Assign significance indicators
  sig = sapply(pvalues.df$p.value, function(y) {

    ifelse(y > 0.01 & y < 0.05, "*",
           ifelse(y > 0.001 & y <= 0.01, "**",
                  ifelse(y <= 0.001, "***", "")
           )
    )

  } )

  pvalues.df = cbind(pvalues.df, sig)

  colnames(pvalues.df)[ncol(pvalues.df)] = ""

  return(pvalues.df)

}

#' sem.model.fits
#' 
#' @keywords internal
#' 
sem.model.fits <- function(...) {

  warning("`sem.model.fits` has been replaced by `rsquared`", call. = FALSE)

}
