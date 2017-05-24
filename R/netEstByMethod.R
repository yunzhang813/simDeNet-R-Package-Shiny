# library(WGCNA)
# library(minet)
# library(bnlearn)


##################
## est.Bayesian ##
##################
est.Bayesian<-function(datExpr){
  genenames <- names(datExpr)
  names(datExpr) <- paste0("V",1:ncol(datExpr))
  ## find linear coeficients from bn.hc
  bn.hc.fit <- bn.fit(hc(datExpr),datExpr)
  coefs <- coef(bn.hc.fit)
  
  ## edge
  child <- parent <- weight <- vector()
  for (i in 1:length(coefs)){
    z <- coefs[[i]][-1]
    child <- c(child, rep(names(coefs)[i],length(z)))
    parent <- c(parent, names(z))
    weight <- c(weight, abs(z))
  }
  child <- gsub("V","C",child)
  parent <- gsub("V","P",parent)
  edge <- data.frame(child, parent, weight)
  
  ## estimated weight on edges 
  nGenes <- ncol(datExpr)
  est.bn0 <- matrix(0, nGenes, nGenes)
  rownames(est.bn0) <- paste0("C",1:nGenes) #row=child
  colnames(est.bn0) <- paste0("P",1:nGenes) #col=parent
  for (i in 1:nrow(edge)){
    est.bn0[child[i],parent[i]] <- weight[i] 
  }
  est.bn <- est.bn0 + t(est.bn0)
  colnames(est.bn) <- rownames(est.bn) <- genenames
  
  return(est.bn)
}


######################
## netEstByMethod() ##
######################

netEstByMethod <- function(data.list,method,softPower=1){
  est.str.list <- vector("list",length=length(data.list))
  names(est.str.list) <- names(data.list)
  if(method=="WGCNA"){
    for(k in 1:length(data.list)){
      est.str.list[[k]] <- adjacency(as.data.frame(t(data.list[[k]])), power=softPower)
    }
  }
  if(method=="ARACNE"){
    for(k in 1:length(data.list)){
      est.str.list[[k]] <- minet(as.data.frame(t(data.list[[k]])), method = "aracne")
    }
  }
  if(method=="Bayesian"){
    for(k in 1:length(data.list)){
      est.str.list[[k]] <- est.Bayesian(as.data.frame(t(data.list[[k]])))
    }
  }
  return(est.str.list)
}
