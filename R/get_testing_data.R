#' Download sample data
#'
#' Helper function to obtain some sample data
#' 
#' @param version Which sample dataset   
#' @param destination Save to this directory
#'
#' @return A list of datasets that can be used to test scGate    
#'
#' @examples
#' \donttest{
#' testing.datasets <- get_testing_data(version = 'hsa.latest')
#' }
#'
#' @export
#'
get_testing_data <- function(version = 'hsa.latest', destination = tempdir()){
  data.folder = file.path(destination,"testing.data")
  if(!dir.exists(data.folder)){
    dir.create(data.folder,recursive = TRUE)
  }
  if(version == 'hsa.latest'){
    testing.data.url = "https://figshare.com/ndownloader/files/31114669?private_link=75b1193bd4c705ffb50b"
    testing.data.path = file.path(data.folder,"testing.dataset.2k.rds")
  }
  if(!file.exists(testing.data.path)){
      download.file(url = testing.data.url,destfile = testing.data.path)
  }
  
  testing.data <- readRDS(testing.data.path)
  return(testing.data)
}
