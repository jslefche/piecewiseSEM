sem.plot = function(
  modelList = NULL, data = NULL, coef.table = NULL, corr.errors = NULL, 
  show.nonsig = TRUE, scaling = 10, alpha = 0.05, ...) {
  
  # Get coefficients
  if(!is.null(modelList) & !is.null(data) & is.null(coef.table))
    
    coef.table = sem.coefs(modelList, data, corr.errors = corr.errors, ...) else

      if(is.null(modelList) & is.null(data) & is.null(coef.table))
        
        stop("Please provide model list and data, or coefficient table!")
      
  # Prepare coef.table
  coef.table$corr.errors = grepl("~~", coef.table[, 1])
  
  coef.table[, 1] = gsub("~~ ", "", coef.table[, 1])
  
  coef.table[, 2] = gsub("~~ ", "", coef.table[, 2])
  
  # Strip transformations
  for(i in 1:2) coef.table[, i] = gsub(".*\\((.*)\\).*", "\\1", coef.table[, i])
  
  # Get vector of labels
  lbls = unlist(coef.table[, 1:2])
  
  lbls = as.character(unname(lbls[!duplicated(lbls)]))
  
  # Shorten label names if necessary
  if(any(sapply(lbls, function(x) nchar(x) > 10))) {
    
    new.lbls = gsub("a|e|i|o|u", "", lbls) } else new.lbls = lbls
  
  if(any(sapply(new.lbls, function(x) nchar(x) > 10))) 
      
      new.lbls = sapply(lbls, function(x) ifelse(nchar(x) > 10, substr(x, 1, 10), x))

  names(new.lbls) = lbls

  # Set graphical parameters
  par(mar = rep(5, 4), xpd = NA)
  
  # Initialize plot
  plot(c(-1.1, 1.1), c(-1.1, 1.1), type = "n", ann = FALSE, axes = FALSE)
  
  # Prepare circle coef.table
  theta = seq(0, 2 * pi, length = 200)
  
  # Add to coef.table.frame
  circle = data.frame(
    x = cos(theta),
    y = sin(theta)
  )
  
  # Get labels
  row.n = seq(1, 200, by = 200/length(lbls))
  
  names(row.n) = lbls
  
  # Get linewidth scale
  if(is.na(scaling)) scl.fctr = rep(1, nrow(coef.table)) else {
    
    scaling = scaling / diff(range(coef.table[, "estimate"]))
    
    scl.fctr = abs(coef.table[, "estimate"]) * scaling
    
  }
  
  # Add arrows
  for(i in 1:nrow(coef.table)) {
    
    resp = row.n[names(row.n) == coef.table[i, 1]]
    
    pred = row.n[names(row.n) == coef.table[i, 2]]
    
    if(show.nonsig == FALSE & coef.table[i, "p.value"] >= alpha) next else {
          lines(x=c(circle[pred, "x"],circle[resp, "x"]),y=c(circle[pred,"y"],circle[resp,"y"])
             , col = ifelse(coef.table[i, "p.value"] <alpha & coef.table[i, "estimate"] > 0, "black", ifelse(coef.table[i, "p.value"] < alpha & coef.table[i, "estimate"] < 0, "red", "grey50"))
             , lwd = scl.fctr[i]
             , lty = ifelse(coef.table[i, "p.value"] < alpha, 1, 2))
    
      arrows(
        x0 = circle[pred, "x"],
        y0 = circle[pred, "y"],
        x1 = circle[resp, "x"]-(circle[resp, "x"]-circle[pred, "x"])/2,
        y1 = circle[resp, "y"]-(circle[resp,"y"]-circle[pred,"y"])/2,
        code = ifelse(coef.table[i, "corr.errors"] == TRUE, 3, 2),
        col = ifelse(coef.table[i, "p.value"] < alpha & coef.table[i, "estimate"] > 0, "black", 
                     ifelse(coef.table[i, "p.value"] < alpha & coef.table[i, "estimate"] < 0, "red", "grey50")),
        lwd = scl.fctr[i],
        lty = ifelse(coef.table[i, "p.value"] < alpha, 1, 2)
      )
      
    }
    
  }
  
  # Plot labels
  for(i in row.n) {
    
    x = circle[i, "x"]
    
    y = circle[i, "y"]
  
    xjust = ifelse(x > 0, 0.5 - (x / 2), 0.5 + (-x / 2))
    
    yjust = ifelse(y > 0, 0, 1)
    
    legend(
      x, 
      y, 
      legend = new.lbls[names(new.lbls) == names(row.n[row.n == i])], 
      cex = 0.8,
      x.intersp = 0,
      xjust = xjust,
      yjust = yjust
      )
    
  }
  
  # Reset par
  par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = FALSE)
  
}
