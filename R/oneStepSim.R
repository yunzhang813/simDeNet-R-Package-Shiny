# library(mvtnorm)
# library(preprocessCore)

##################
## oneStepSim() ##
##################

oneStepSim <- function(n.samp, mu.T, mu.N, Sigma.T=NULL, Sigma.N=NULL, prop.T="seq",
                       # structure for Sigma.T
                       block.size, rho, dd=NULL, str.type="interchangeable",
                       # selected genes to add structure
                       select.gene="first"){
  ## parameters
  m.gene <- length(mu.T)
  if(is.null(Sigma.N)){Sigma.N <- diag(m.gene)} #identity matrix for Sigma.N
  if(is.null(Sigma.T)){
    out.Sigma.T <- generate.Sigma(m.gene=m.gene, dd=dd, rho=rho, block.size=block.size, str.type=str.type)
    Sigma.T <- out.Sigma.T$Sigma
    true.str.T <- out.Sigma.T$true.str
  }
  
  ## selected genes to add structure
  if(select.gene=="first"){oo <- 1:m.gene}
  if(select.gene=="random"){oo <- sample(1:m.gene, m.gene, replace=FALSE)}
  if(select.gene=="DEG"){oo <- order(abs(mu.T-mu.N), decreasing = TRUE)}
  if(select.gene=="non-DEG"){oo <- order(abs(mu.N-mu.T), decreasing = FALSE)}
  
  ## re-order genes
  mu.T <- mu.T[oo]
  mu.N <- mu.N[oo]
  
  ## simulate expression for pure and mixed samples
  out.expr <- sim.expr(n.samp=n.samp, mu.T=mu.T, mu.N=mu.N, Sigma.T=Sigma.T, Sigma.N=Sigma.N, prop.T=prop.T)
  
  return(list("m.gene"=m.gene, "Sigma.T"=Sigma.T, "Sigma.N"=Sigma.N, "true.str.T"=true.str.T,
              "expr.pure.T"=out.expr$expr.pure.T, "expr.pure.N"=out.expr$expr.pure.N, 
              "expr.mixed"=out.expr$expr.mixed, "true.prop"=out.expr$true.prop))
}


################
## sim.expr() ##
################

sim.expr <- function(n.samp, mu.T, mu.N, Sigma.T, Sigma.N, prop.T="seq"){
  ## check
  if(length(mu.T)!=length(mu.N)){stop("Error: mu.T and mu.N should be have the same length.")}
  if(sum(dim(Sigma.T)==length(mu.T))!=2){stop("Error: Dimension of Sigma.T should match length of mu.T.")}
  if(sum(dim(Sigma.N)==length(mu.N))!=2){stop("Error: Dimension of Sigma.N should match length of mu.N.")}
  if((prop.T!="seq")&&(prop.T!="unif")&(length(prop.T)!=n.samp)){stop("Error: customized prop.T should have length of n.samp.")}
  if((prop.T!="seq")&&(prop.T!="unif")&((min(prop.T)<0)|(max(prop.T)>1))){stop("Error: customized prop.T should range from 0 to 1.")}
  
  ## pure samples
  expr.pure.T <- t(rmvnorm(n.samp, mu.T, Sigma.T))
  expr.pure.N <- t(rmvnorm(n.samp, mu.N, Sigma.N))
  expr.pure.T[expr.pure.T<0] <- 0
  expr.pure.N[expr.pure.N<0] <- 0
  
  ## set prop.T
  if(is.character(prop.T)&&prop.T=="seq"){prop.T <- seq(0, 1, length=n.samp)}
  else if(is.character(prop.T)&&prop.T=="unif"){prop.T <- runif(n.samp, 0, 1)}
  
  ## mixed samples
  design.mixed <- rbind(diag(prop.T),diag(1-prop.T))
  ## mixing in raw scale
  expr.pure.raw <- 2^cbind(expr.pure.T, expr.pure.N)
  expr.mixed.raw <- expr.pure.raw %*% design.mixed
  ## mixed expression
  expr.mixed <- log2(expr.mixed.raw)
  
  ## quantile normalization
  expr.pure.T <- normalize.quantiles(expr.pure.T)
  expr.pure.N <- normalize.quantiles(expr.pure.N)
  expr.mixed <- normalize.quantiles(expr.mixed)
  
  ## sample names in columns
  colnames(expr.pure.T) <- paste0("pure.T.samp",1:n.samp)
  colnames(expr.pure.N) <- paste0("pure.N.samp",1:n.samp)
  colnames(expr.mixed) <- paste0("mixed.samp",1:n.samp)
  ## gene names in rows
  if(!is.null(names(mu.T))){
    rownames(expr.pure.T) <- rownames(expr.pure.N) <- rownames(expr.mixed) <- names(mu.T)
  }
  
  return(list("expr.pure.T"=expr.pure.T, "expr.pure.N"=expr.pure.N, "expr.mixed"=expr.mixed, "true.prop"=prop.T))
}


####################
## generate.Sigma ##
####################

generate.Sigma <- function(m.gene, block.size, rho, dd=NULL, str.type="interchangeable"){
  ## check
  if(sum(block.size)>m.gene){stop("Error: sum of block.size must be smaller than m.gene.")}
  if(length(rho)!=length(block.size)){stop("Error: rho and block.size must be have the same length.")}
  if(!is.null(dd)&length(dd)!=m.gene){stop("Error: dd must have length equal to m.gene.")}
  if(length(str.type)!=length(block.size)&length(str.type)!=1){stop("Error: either specify one str.type for all the blocks or specify different str.type for each block.")}
  
  ## weight matrix
  if(is.null(dd)){dd <- rep(1,m.gene)}
  Wmat <- diag(dd)
  ##Sigma and truth matrix
  Sigma <- truth <- matrix(0, m.gene, m.gene)
  
  ## block structure
  if (length(str.type)==1){str.type <- rep(str.type,length(rho))}
  block.ends <- cumsum(block.size)
  block.starts <- block.ends-block.size+1
  
  for (k in 1:length(block.size)){
    ## interchangeable correlation
    if(str.type[k]=="interchangeable"){
      block.k <- matrix(rho[k], block.size[k], block.size[k])
      diag(block.k) <- 0
      ## add structure to Sigma
      Sigma[block.starts[k]:block.ends[k],block.starts[k]:block.ends[k]] <- Sigma[block.starts[k]:block.ends[k],block.starts[k]:block.ends[k]] + block.k
      ## truth
      if(rho[k]!=0){truth[block.starts[k]:block.ends[k],block.starts[k]:block.ends[k]] <- truth[block.starts[k]:block.ends[k],block.starts[k]:block.ends[k]] + matrix(as.numeric(block.k==rho[k]), block.size[k], block.size[k])}
    }
    
    ## decaying correlation
    if(str.type[k]=="decaying"){
      block.k <- matrix(NA, block.size[k], block.size[k])
      for(i in 1:block.size[k]){
        for(j in 1:block.size[k]){
          block.k[i,j] <- rho[k]^abs(i-j) 
        }
      }
      diag(block.k) <- 0
      ## add structure to Sigma
      Sigma[block.starts[k]:block.ends[k],block.starts[k]:block.ends[k]] <- Sigma[block.starts[k]:block.ends[k],block.starts[k]:block.ends[k]] + block.k
      ## truth
      if(rho[k]!=0){truth[block.starts[k]:block.ends[k],block.starts[k]:block.ends[k]] <- truth[block.starts[k]:block.ends[k],block.starts[k]:block.ends[k]] + matrix(as.numeric(block.k==rho[k]), block.size[k], block.size[k])}
    }
    
    ## star-shape correlation
    if(str.type[k]=="star"){
      block.k <- matrix(rho[k], block.size[k], block.size[k])
      for(i in 2:block.size[k]){
        for(j in 2:block.size[k]){
          block.k[i,j] <- rho[k]^2 
        }
      }
      diag(block.k) <- 0
      ## add structure to Sigma
      Sigma[block.starts[k]:block.ends[k],block.starts[k]:block.ends[k]] <- Sigma[block.starts[k]:block.ends[k],block.starts[k]:block.ends[k]] + block.k
      ## truth
      if(rho[k]!=0){truth[block.starts[k]:block.ends[k],block.starts[k]:block.ends[k]] <- truth[block.starts[k]:block.ends[k],block.starts[k]:block.ends[k]] + matrix(as.numeric(block.k==rho[k]), block.size[k], block.size[k])}
    }
  }
  ## more cares on the diagnal
  diag(Sigma) <- 1
  diag(truth) <- 0
  
  ## Wmat %*% Sigma %*% Wmat
  Sigma.final <- Wmat %*% Sigma %*% Wmat
  
  return(list("Sigma"=Sigma.final, "true.str"=truth))
}




