model.to.table <- function(scGate.model){
  tab <- data.frame(levels=NULL,use_as = NULL,name = NULL,signature = NULL)
  
  ### Define first column: Levels
  levels <- scGate.model%>%names()
  lev <- rep(levels,rep(2,length(levels)))
  
  len_pos_neg <- lapply(scGate.model,function(x){ 
    res = lapply(x,function(y){
      length(y)
    })
    return(res)
  })%>%unlist()
  
  extended.levels <- rep(lev,len_pos_neg)
  
  # second column: "use_as"
  useas <- rep(c("positive","negative"),length(levels))
  useas <- rep(useas , len_pos_neg)
  
  #third column: name
  signature.names <- lapply(scGate.model,function(x){ 
    res = lapply(x,function(y){
      names(y)
    })
    return(res)
  })%>%unlist()
  
  ## Four column: signature
  signatures <- lapply(scGate.model,function(x){ 
    res = lapply(x,function(y){
      lapply(y, function(z){
        paste(z,collapse = ";")
      })
    })
    return(res)
  })%>%unlist()
  
  tab <- data.frame("levels"=extended.levels,"use_as" = useas, "name" = signature.names,"signature" = signatures)
  return(tab)
}
