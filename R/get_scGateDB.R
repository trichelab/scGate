
#' Load scGate model database
#'
#' Download, update or load local version of the scGate model database.
#' These are stored in a GitHub repository, whence you can download specific 
#' versions of the database.
#' 
#' @param destination Destination path for storing the DB. Default is tempdir(); if you wish to edit locally the models and
#'    link them to the current project, set this parameter to a new directory name, e.g. scGateDB
#' @param force_update  Whether to update an existing database.
#' @param version Specify the version of the scGate_models database (e.g. 'v0.1'). By default downloads the latest available version.
#' @param repo_url  URL path to scGate model repository database
#' @param branch  branch of the scGate model repository, either 'master' (default) or 'dev' for the latest models 
#' @param verbose  display progress messages
#'
#' @return A list of models, organized according to the folder structure of the database. See the examples below.
#'
#' @details
#' Models for scGate are dataframes where each line is a signature
#' for a given filtering level. A database of models can be downloaded
#' using the function \code{get_scGateDB}. You may directly use the models from
#' the database, or edit one of these models to generate your own custom gating
#' model.  
#'
#' @examples
#' scGate.model.db <- get_scGateDB()
#' # To see a specific model, browse the list of models:
#' scGate.model.db$human$generic$Myeloid
#' # Apply scGate with this model
#' data(query.seurat)
#' query <- scGate(query.seurat, model=scGate.model.db$human$generic$Myeloid, reduction="pca")
#'
#' @seealso \code{\link{scGate}} \code{\link{load_scGate_model}}
#'
#' @importFrom dplyr %>%  
#' @importFrom utils download.file unzip read.table
#'
#' @export
#'
get_scGateDB <- function(destination = tempdir(),
                         force_update = FALSE,
                         version = "latest",
                         branch=c("master","dev"), 
                         verbose=FALSE,
                         repo_url="https://github.com/carmonalab/scGate_models"){
  
  branch = branch[1]
  if (version == "latest") {
    repo_url_zip = sprintf("%s/archive/%s.zip", repo_url,branch)
    repo.name <- paste0("scGate_models-",branch)
    repo.name.v <- repo.name
  } else {
    repo_url_zip = sprintf("%s/archive/refs/tags/%s.zip", repo_url, version)
    repo.name = sprintf("scGate_models-%s", version)
    #for some reason GitHub remove the 'v' from repo name after unzipping
    repo.name.v <- sprintf("scGate_models-%s", gsub("^v","",version, perl=TRUE)) 
  }
  destination <- normalizePath(destination, winslash = "/")
  repo_path = file.path(destination,repo.name)
  repo_path.v = file.path(destination,repo.name.v)
  temp <- tempfile()
  
  if(!dir.exists(repo_path)){
    if(!dir.exists(destination)) {
      dir.create(destination)
    }
    download.file(repo_url_zip,temp)
    unzip(temp,exdir = destination)
    unlink(temp)
  }else if(force_update){
    download.file(repo_url_zip,temp)
    system(sprintf("rm -r %s",repo_path))  # this ensure that db would be completely overwritten and old model will not persist. 
    unzip(temp,exdir = destination, overwrite = force_update)
    unlink(temp)
  }else{
    message(sprintf("Using local version of repo %s. If you want update it, set option force_update = TRUE",repo.name))
  }
  
  #Now load the models into a list structure
  allfiles <- list.files(repo_path.v, recursive = TRUE)
  modelfiles <- grep("scGate_Model.tsv", allfiles, value = TRUE)
  uniq_dirs <- sort(unique(dirname(modelfiles)))
  
  model_db <- list()
  for (dir in uniq_dirs) {
    sub <- strsplit(dir, split="/")[[1]]
    model_path <- file.path(repo_path.v, dir)
    
    if (length(sub)==0) {
      stop("Error in scGate DB format")
    } else if (length(sub)==1) {
      if(verbose) message(paste("loading ",model_path))
      model_db[[sub[1]]] <- load.model.helper(model_path,verbose=verbose)
    } else if (length(sub)==2) {
      if(verbose) message(paste("loading ",model_path))
      model_db[[sub[1]]][[sub[2]]] <- load.model.helper(model_path,verbose=verbose)
    } else if (length(sub)==2) {
      if(verbose) message(paste("loading ",model_path))
      model_db[[sub[1]]][[sub[2]]][[sub[[3]]]] <- load.model.helper(model_path,verbose=verbose)
    } else {
      message(sprintf("Warning: max depth for scGate models is 3. Skipping folder %s", model_path))
    } 
  }
  return(model_db)
}
