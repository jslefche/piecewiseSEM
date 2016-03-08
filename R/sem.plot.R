sem.plot = function(modelList = NULL, data = NULL, table = NULL, ...) {
  
  # Get coefficients
  if(!is.null(modelList) & !is.null(data) & is.null(table))
    
    table = sem.coefs(modelList, data, ...) else

      if(is.null(modelList) & is.null(data) & is.null(table))
        
        stop("Please provide model list and data, or coefficient table!")
      
  
  # Get vector of labels
  lbls = unlist(table[, 1:2])
  
  lbls = as.character(unname(lbls[!duplicated(lbls)]))

  # Set graphical parameters
  par(mar = rep(3, 4), xpd = NA)
  
  # Initialize plot
  plot(c(-1.1, 1.1), c(-1.1, 1.1), type = "n", ann = FALSE, axes = FALSE)
  
  # Prepare circle table
  theta = seq(0, 2 * pi, length = 200)
  
  # Add to table.frame
  circle = data.frame(
    x = cos(theta),
    y = sin(theta)
  )
  
  # Get labels
  row.n = seq(1, 200, by = 200/length(lbls))
  
  names(row.n) = lbls
  
  # Get linewidth scale
  scaling = 10 / diff(range(table[, "estimate"]))
  
  scl.fctr = abs(table[, "estimate"]) * scaling
  
  # Add arrows
  for(i in 1:nrow(table)) {
    
    resp = row.n[names(row.n) == table[i, 1]]
    
    pred = row.n[names(row.n) == table[i, 2]]
    
    arrows(
      x0 = circle[pred, "x"],
      y0 = circle[pred, "y"],
      x1 = circle[resp, "x"],
      y1 = circle[resp, "y"],
      col = ifelse(table[i, "p.value"] < 0.05 & table[i, "estimate"] > 0, "black", 
                   ifelse(table[i, "p.value"] < 0.05 & table[i, "estimate"] < 0, "red", "grey50")),
      lwd = scl.fctr[i],
      lty = ifelse(table[i, "p.value"] < 0.05, 1, 2)
    )
    
  }
  
  # Plot labels
  for(i in row.n) {
    
    x = circle[i, "x"]
    
    y = circle[i, "y"]
  
    xjust = ifelse(x > 0, -0.3, 1.3)
    
    yjust = ifelse(y > 0 & x > 0, 0, 1)
    
    legend(
      x, 
      y, 
      legend = names(row.n[row.n == i]), 
      cex = 0.8,
      x.intersp = 0,
      xjust = xjust,
      yjust = yjust
      )
    
  }
  
}