#' Performance metrics
#'
#' Evaluate model performance for binary tasks
#' 
#' @param actual Logical or numeric binary vector giving the actual cell labels.
#' @param pred  Logical or numeric binary vector giving predicted cell labels. 
#' @param return_contingency  Return the contingency table? (FALSE) 
#'
#' @return Performance metrics (Precision, Recall, MCC) for actual vs predicted
#'
#' @examples
#' results <- performance.metrics(actual= sample(c(1,0),20,replace=TRUE),
#'     pred =  sample(c(1,0),20,replace=TRUE,prob = c(0.65,0.35) ) )
#'
#' @export
#'
performance.metrics <- function(actual, pred, return_contingency=FALSE){

  actual <- as.numeric(actual +0)  
  pred <- as.numeric(pred +0)  
  tp <- sum(actual&pred)
  tn <- sum((!actual)&(!pred))
  fn <- sum(actual&(!pred))
  fp <- sum((!actual)&pred)  
  
  PREC <- tp/(tp +fp)
  REC <- tp/(tp + fn)
  #sqrt_ <- sqrt((tp + fp)*(tp+fn)*(tn+fp)*(tn+fn))
  sqrt_ <- exp(0.5* sum(log(c(tp+fp, tp+fn, tn+fp, tn+fn))) )
  MCC <- (tp*tn - fp*fn) / sqrt_
  
  
  if (return_contingency) {  
    
    ct <- table(actual,pred)
    ## ordering contingency table, 
    ## but avoiding errors when all predictions (or all actual cells) are equals
    nam.act <- as.character(sort(unique(actual), decreasing = TRUE))
    nam.pred <- as.character(sort(unique(pred), decreasing = TRUE))
    ct <- ct[nam.act,nam.pred]
    res <- list('counting' = ct, 'summary' = res.Summary )
    return(res)
    
  } else {

    res.Summary <- c(PREC,REC,MCC)
    names(res.Summary) <- c("PREC","REC","MCC")
    return(res.Summary)
  
  }
  
}
