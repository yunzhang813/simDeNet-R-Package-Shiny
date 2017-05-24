## library(ISOpureR)

############
## deconv ##
############

deconv <- function(mixed, ref, seed=123, ...){
  set.seed(seed)
  ## transform the expression data to raw scale
  mixed.dat <- 2^mixed
  ref.dat <- 2^ref
  
  ## step 1
  ISOpureS1model <- ISOpure.step1.CPE(mixed.dat, ref.dat, ...)
  ## step 2
  ISOpureS2model <- ISOpure.step2.PPE(mixed.dat, ref.dat, ISOpureS1model)
  dat.ISOpure <- ISOpureS2model$cc_cancerprofiles # decovoluted non-ref samples (at the raw scale)
  est.prop <- ISOpureS2model$alphapurities # estimated proportion of non-ref samples
  
  ## output log2-transformed expression
  return(list("expr.deconv"=log2(dat.ISOpure), "est.prop"=est.prop))
}
