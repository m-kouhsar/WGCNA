
moduleEigengenes_ <- function(expr,colors){
  m = unique(colors)
  me = as.data.frame(matrix(data = NA , ncol = length(m) , nrow = ncol(expr)))
  names(me) = m
  rownames(me) = colnames(expr)
  
  for (i in m) {
    m.expr = expr[colors==i , ]
    pca = prcomp(t(m.expr),center = T , scale. = T)
    me[,i] = pca$x[,1]
  }
  names(me)=paste0("ME",names(me))
  return(me)
}