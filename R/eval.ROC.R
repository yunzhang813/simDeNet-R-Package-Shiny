# library(ROCR)

##############
## eval.ROC ##
##############


eval.ROC <- function(est.str, true.str, plot.ROC=TRUE, show.AUC=TRUE, zoom=TRUE, lightPDF=FALSE, ...){
  # if(sum(dim(est.str)!=dim(est.str))>0){stop("Error: est.str and true.str must have the same dimemsion.")}
  # if(!isSymmetric(unname(est.str))|!isSymmetric(unname(true.str))){stop("Error: est.str and true.str must be symmetric.")}
  
  ## evaluate ROC and AUC
  ind <- lower.tri(est.str)
  pred <- prediction(est.str[ind], true.str[ind])
  perf <- performance(pred, "tpr", "fpr")
  auc <- performance(pred, "auc")
  auc.value <- auc@y.values[[1]]
  
  ## parameters for re-size ROC plot
  n.unique <- length(unique(est.str[ind]))
  xval <- slot(perf,"x.values")[[1]][n.unique-1]
  yval <- slot(perf,"y.values")[[1]][n.unique-1]
  
  ## plot
  if(plot.ROC){
    perf.plot <- perf
    if(lightPDF){
      set.seed = 123
      sparsePoints <- sort(sample(n.unique, 10000))
      slot(perf.plot,"x.values")[[1]] <- slot(perf.plot,"x.values")[[1]][sparsePoints]
      slot(perf.plot,"y.values")[[1]] <- slot(perf.plot,"y.values")[[1]][sparsePoints]
    }
    if(zoom){plot(perf.plot, xlim=c(0,xval), ylim=c(0,yval), ...)}
    else plot(perf.plot, ...)
    if(show.AUC){legend("bottom",paste0("AUC=",round(auc.value,3)),bty = "n")}
  }
  
  return(list("perf"=perf, "AUC"=auc.value, "xval"=xval, "yval"=yval))
}

