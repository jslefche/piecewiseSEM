#' Plotting of Piecewise Structural Equation Models
#' 
#' @description plot.psem uses [DiagrammeR] to generate path diagrams 
#' of `piecewiseSEM`` fits within R.
#' 
#' @param x a [psem()] object
#' @param return whether to return the output from [DiagrammeR::create_graph()] for modification and later plotting
#' @param node_attrs List of node attributes to override defaults of rectangular nodes with black outline and white fill. See [here](http://visualizers.co/diagrammer/articles/node-edge-data-frames.html) and [here](http://visualizers.co/diagrammer/articles/graphviz-mermaid.html) for a more complete rundown of options.
#' @param edge_attrs List of edge attributes to override defaults of solid black arrows. See [here](http://visualizers.co/diagrammer/articles/node-edge-data-frames.html) and [here](http://visualizers.co/diagrammer/articles/graphviz-mermaid.html) for a more complete rundown of options.
#' @param ns_dashed If TRUE, paths that are not different from 0 will be dashed rather than solid, unless the whole is overridden in `edge_attrs`
#' @param alpha The alpha level for assessing whether a path is different from 0
#' @param show What types of path coefficients are shown? Default `"std"` is standardized coefficients. For undstandardized, use `"unstd"`
#' @param digits How many significant digits should be shown?
#' @param add_edge_label_spaces Should spaces by added on either side of edge labels? Default is `TRUE` as otherwise paths too often overlap edges.
#' @param ... Other arguments to [DiagrammeR::render_graph()]
#' 
#' @return Returns an object of class [DiagrammeR::dgr_graph]
#' 
#' @author Jarrett Byrnes <jarrett.byrnes@@umb.edu>
#' 
#' @examples 
#' data(keeley)
#' 
#' mod <- psem(
#'   lm(rich ~ cover, data=keeley),
#'   lm(cover ~ firesev, data=keeley),
#'   lm(firesev ~ age, data=keeley),
#'   data = keeley
#' )
#' 
#' plot(mod)
#' 
#' ### More customized plot
#' 
#' plot(mod, node_attrs = list(
#'   shape = "rectangle", color = "black",
#'   fillcolor = "orange", x = 3, y=1:4))
#'   
#' @import DiagrammeR
#' 
#' @method plot psem
#' 
#' @export
#' 
plot.psem <- function(x, return=FALSE,
                      node_attrs = data.frame(shape = "rectangle", color = "black",
                                              fillcolor = "white"),
                      edge_attrs = data.frame(style = "solid", color="black"),
                      ns_dashed = T, alpha=0.05,
                      show = "std", digits = 3, 
                      add_edge_label_spaces = TRUE, ...
                      ){
  
  #get the coefficients table
  ctab <- coefs(x)
  ctab$Response <- as.character(ctab$Response)
  ctab$Predictor <- as.character(ctab$Predictor)
  
  #make a nodes DF
  unique_nodes <- unique(c(ctab$Response, ctab$Predictor))
  nodes <- create_node_df(n = length(unique_nodes),
                          nodes = unique_nodes,
                          type = "lower",
                          label = unique_nodes)
  nodes <- cbind(nodes, node_attrs)
  nodes[] <- lapply(nodes, as.character)
  nodes$id <- as.numeric(nodes$id)

  #make an edges DF
  edges <- create_edge_df(
                          from = match(ctab$Predictor, unique_nodes),
                          to = match(ctab$Response, unique_nodes))
  
  edges <- data.frame(edges, edge_attrs)
  edges[] <- lapply(edges, as.character)
  edges$id <- as.numeric(edges$id)
  edges$from <- as.numeric(edges$from)
  edges$to <- as.numeric(edges$to)
  if(ns_dashed) edges$style[which(ctab$P.Value>alpha)] <- "dashed"
  if(show == "std") edges$label = round(ctab$`Std.Estimate`, digits)
  if(show == "unstd") edges$label = round(ctab$Estimate, digits)
  if(add_edge_label_spaces) edges$label = paste0(" ", edges$label, " ")
  
  #turn into a graph
  sem_graph <- create_graph(nodes, edges, directed=TRUE)
  
  if(return) return(sem_graph)
  
  render_graph(sem_graph, ...)

}