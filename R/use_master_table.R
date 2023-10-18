## This function allows to complete signatures in a table based model by using the name signature and a provided master.table of signatures
# the master.table must be a two column data.frame with two columns : 1) name: contains the signature names and 
# 2)signature: this column contain the genes present in each signature (separated with a semicolon) 
use_master_table <- function(df.model, master.table, name = "name",descript = "signature"){
  
  ## First, check if descript field exists in input table
  if(!descript %in% colnames(df.model)){
    df.model[descript] <- ""
  }
  
  ## second, check if there is something to do (or exit)
  input.sign <- df.model[[descript]]
  input.names <- df.model[[name]]
  complete.from.master.table <- as.vector(is.null(input.sign)|input.sign == "" | is.na(input.sign))
  
  if(!any(complete.from.master.table)){
    return(df.model)
  }
  
  #Does the master table exist?
  if(!file.exists(master.table)){
    stop("master_table.tsv file must be present in your 'model folder' unless signatures are completely specified in the models")
  }
  master.table <- read.table(master.table, sep ="\t", header =T) 
  
  # sanity check:
  warn <- setdiff(input.names[complete.from.master.table],master.table[[name]])
  if(length(warn)>0){
    stop(sprintf("signatures '%s' are not present in the provided master.signature table",paste(warn,collapse = " ; " )))
  }
  
  merged <- merge(df.model[complete.from.master.table,],master.table,by = name,suffixes = c("","_from.master"),all.x = T)
  vect <- merged[[paste0(descript,"_from.master")]]
  names(vect) <- merged[[name]]
  if(merged%>%nrow > sum(complete.from.master.table)){
    stop("please, master.table must not contain duplicated signature names")
  }
  
  # replace ensuring correct order (notice that the merge can change row order)
  nms <- df.model[complete.from.master.table,name]
  df.model[complete.from.master.table,descript] <- vect[nms]
  #  df.model[descript] <- output.sign
  return(df.model)
}
