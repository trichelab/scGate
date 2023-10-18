load.model.helper <- function(models_path, master.table = "master_table.tsv",  verbose=verbose) {

  df.models.toimpute <- list()
  files.to.impute <- list.files(file.path(models_path),"_scGate_Model.tsv")
  if(length(files.to.impute)==0){
    stop("Please, provide some model table files in your 'model folder' or set models_path = NULL for using the default ones")
  }
  # load models to impute
  for(f in files.to.impute){
    model.name <- strsplit(f,"_scGate_Model.tsv")[[1]][1]
    if(verbose) message(paste0("loading ",model.name))
    df.models.toimpute[[model.name]] <- read.table(file.path(models_path,f),sep ="\t",header =T)
  }
  # signature imputing
  imputed.models <-  lapply(df.models.toimpute,function(df){
    use_master_table(df.model = df, master.table = file.path(models_path, master.table))
  })
  model.list <- imputed.models
      

  return(model.list)
}
