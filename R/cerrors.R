cerror <- function(.formula, modelList, data) {

  if(class(.formula) == "formula.cerror") {

    x <- pcor(.formula, modelList, data)

    p <- pcortest(x)

  }


  return(x, p)

}

pcor <- function() {

  if(class(.formula) == "formula.cerror")

    vars <- gsub(" " , "", unlist(strsplit(.formula, "~~"))) else

      vars <- gsub(" ", "", unlist(strsplit(deparse(.formula), "~")))



  vars <- unlist(strsplit(deparse(.formula), "~"))

}

pcortest <- function() {

}


`%~~%` <- function(e1, e2) {

  x <- paste(deparse(substitute(e1)), "~~", deparse(substitute(e2)))

  # x <- call(x)

  class(x) <- "formula.cerror"

  return(x)

}
