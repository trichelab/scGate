#' Test your model
#'
#' Wrapper for fast model testing on 3 sampled datasets 
#' 
#' @param model scGate model in data.frame format 
#' @param testing.version  Character indicating the version of testing tatasets
#'     to be used. By default "hsa-latest" will be used. It will be ignored if
#'     a custom dataset is provided (in Seurat format). 
#' @param custom.dataset  Seurat object to be used as a testing dataset. For
#'     testing purposes, metadata seurat object must contain a column named
#'     'cell_type' to be used as a gold standard. Also a set of positive
#'     targets must be provided in the target variable. 
#' @param target Positive target cell types. If default testing version is used
#'     this variable must be a character indicating one of the available target
#'     models ('immune','Lymphoid','Myeloid','Tcell','Bcell','CD8T','CD4T',
#'     'NK','MoMacDC','Plasma_cell','PanBcell'). 
#'     If a custom dataset is provided in Seurat format, this variable must be
#'     a vector of positive cell types in your data. The last case also require
#'     that such labels were named as in your cell_type meta.data column. 
#' @param plot Whether to return plots to device
#' @return Returns performance metrics for the benchmarking datasets, and optionally
#'     plots of the predicted cell type labels in reduced dimensionality space. 
#' @examples
#' \donttest{
#' scGate.model.db <- get_scGateDB()
#' # Browse the list of models and select one:
#' model.panBcell <-  scGate.model.db$human$generic$PanBcell
#' # Test the model with available testing datasets
#' panBcell.performance <- test_my_model(model.panBcell, target = "PanBcell")
#' model.Myeloid <-  scGate.model.db$human$generic$Myeloid
#' myeloid.performance <- test_my_model(model.Myeloid, target = "Myeloid")
#' }     
#' @importFrom utils download.file
#' @importFrom methods is
#' @importFrom patchwork wrap_plots
#' @export

test_my_model <- function(model, testing.version = 'hsa.latest',
                          custom.dataset = NULL,target = NULL,
                          plot = TRUE){
  
  performance.computation  <- ifelse (is.null(target), FALSE, TRUE)
  
  if (is(custom.dataset, "Seurat")){
    testing.datasets <- list()
    testing.datasets[["custom.dataset"]] <- custom.dataset
    custom <- TRUE
  } else { 
    custom <- FALSE
  }

  if(!custom){
    targets <- c('immune','Lymphoid','Myeloid','Tcell','Bcell','CD8T','CD4T','NK','MoMacDC','Plasma_cell','PanBcell')
    
    if(is.null(target)){
      message("warning: target cell_type not provided. Avoiding performance computation")  
      performance.computation <- FALSE
    }else if(!target %in% targets){
      stop(sprintf("target must be one of %s; or NULL for avoiding performance computation",paste(targets,collapse = "';'")))
    }
    
    ## check dataset version
    available.datasets = c("hsa.latest")
    if(!testing.version %in% available.datasets){
      stop("Please provide a valid testing.version parameter or provide a custom.dataset in seurat format")
    }
    
    # load testing datasets
    if(testing.version == "hsa.latest"){
      testing.datasets <- get_testing_data(version = testing.version)
    }
  }  
  
  if(custom){
    if(!"cell_type" %in% colnames(custom.dataset@meta.data)){
      stop("please, provide a 'cell_type' column to be used as reference cell type")
    }
    
    if(is.null(target)){
      message("warning: target cell_type not provided. Avoiding performance computation")  
      performance.computation <- FALSE
    }else if(any(!target %in% custom.dataset$cell_type)){
      stop("all target celltypes must be included in cell_type metadata field. Otherwise, set target = NULL for avoiding performance computation")
    }
  }
  
  plt.out <- list()
  perf.out <- list()
  output <- list()
  # Drop is.pure cols if exists
  for(dset in names(testing.datasets)){
    obj <- testing.datasets[[dset]]
    plt <- list()
    cols <- colnames(obj@meta.data)
    dropcols = grep("^is.pure",cols,value =TRUE) %>% unique()
    if(length(dropcols)>0){
      for(col in dropcols){
        obj[[col]] <- NULL   
      }
    }
    
    ## scGate filtering
    obj <- scGate(obj, model = model, assay = DefaultAssay(obj))

    # add annotation plot
    nname <- sprintf("%s manual annot",dset)
    plt <- DimPlot(obj, group.by = "cell_type", label = TRUE,
                   repel =TRUE, label.size = 3) + 
      ggtitle(nname) + NoLegend() +  theme(aspect.ratio = 1)

    # add one DimPlot by model level
    pure.plot <- DimPlot(obj, group.by = "is.pure", cols = list("Pure"="green","Impure"="gray")) +
      theme(aspect.ratio = 1)
    plt <- list("Annotation"=plt, "Gating"=pure.plot)
    
    #reserve plots of this dset
    plt.out[[dset]] <- patchwork::wrap_plots(plt,ncol = length(plt))
    
    if(performance.computation){
      if(!custom){    
        performance = scGate::performance.metrics(actual = obj@meta.data[,target],
                                                  pred = obj$`is.pure`== "Pure")
      }else{
        performance = scGate::performance.metrics(actual = obj@meta.data$cell_type %in% target,
                                                  pred = obj$`is.pure`== "Pure")
      }
      perf.out[[dset]] <- performance 
    }
    output[[dset]] <- obj
    
  }
  
  if(performance.computation){
    
    if(!custom){
      perf <- Reduce(rbind,perf.out)
      rownames(perf) <- names(perf.out)  
    } else {
      perf <- perf.out
    }
  }
  
  if(plot) {
    print(patchwork::wrap_plots(plt.out, ncol = 1))
  }

  if(performance.computation){
    return(list(performance = perf, plots = plt.out, objects = output))
  }else{
    return(list(plots = plt.out, objects = output))
  }
}
