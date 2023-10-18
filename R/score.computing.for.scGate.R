score.computing.for.scGate <- function(data, model, bpp=SerialParam(), assay="RNA", slot="data",
                                       add.sign=NULL, keep.ranks=FALSE, maxRank=1500) {
  
  comb <- dplyr::bind_rows(model, .id = "Model_ID")
  # extract unique signatures
  model.uniq <- comb %>% dplyr::distinct(.data$name, .data$signature, .keep_all = T) 
  
  #Stop if there are signature with same name but different genes
  t <- table(model.uniq$name)
  dup <- names(t[t>1])
  if (length(dup)>0) {
    s <- paste(dup,collapse=", ")
    stop(sprintf("Different gene sets have been assigned to signature with same name: %s", s))
  }
  
  ## generate list object to be used in computing stage
  signatures <- model.uniq$signature %>% strsplit("[,; ]+") %>% lapply(unlist)   #also allow comma or space
  names(signatures) <- model.uniq$name
  
  if (!is.null(add.sign)) {
    if (!inherits(add.sign, "list")) {
      add.sign <- list("Additional_signature"=add.sign)
    }
    signatures <- append(signatures, add.sign)
  }
  
  data <- AddModuleScore_UCell(data, features = signatures, assay=assay, slot=slot,
                                      BPPARAM = bpp, storeRanks = keep.ranks, maxRank = maxRank)
  
  return(data)
}
