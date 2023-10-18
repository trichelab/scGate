#' Load a single scGate model
#'
#' Loads a custom scGate model into R. For the format of these models, have a
#' view/edit one of the default models obtained with \code{\link{get_scGateDB}}
#'
#' @param model_file scGate model file, in .tsv format.
#' @param master.table File name of the master table (in repo_path folder) that contains cell type signatures.
#'
#' @return A scGate model in dataframe format, which can given as input to the \code{\link{scGate}} function.
#'
#' @examples
#' dir <- tempdir() # this may also be set to your working directory
#' models <- get_scGateDB(destination=dir)
#' # Original or edited model
#' model.path <- paste0(dir,"/scGate_models-master/human/generic/Bcell_scGate_Model.tsv")
#' master.path <- paste0(dir,"/scGate_models-master/human/generic/master_table.tsv")
#' my.model <- load_scGate_model(model.path, master.path)
#' my.model
#'
#' @seealso \code{\link{scGate}} \code{\link{get_scGateDB}} 
#'
#' @importFrom utils read.table
#'
#' @export
#'
load_scGate_model <- function(model_file, master.table = "master_table.tsv") {
  
  model <- read.table(model_file, sep ="\t",header =TRUE)
  model <- use_master_table(model, master.table = master.table)
  
  return(model)

}
