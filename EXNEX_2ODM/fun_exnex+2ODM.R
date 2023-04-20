
exnex.2ODM <- function(J           = 6,
                       n           = rep(30, J),
                       n.exch      = 1,
                       nex.mean    = -1,
                       nex.prec    = 1/3,
                       ex.mu.mean  = c(1),
                       ex.mu.prec  = c(1/3),
                       ex.tau.prec = c(1),
                       dist.prior  = NULL,
                       p0          = c(0.15, 0.2, 0.3, 0.4, 0.5, 0.6),
                       H0          = rep(0.2, J),
                       H1          = rep(0.4, J),
                       tau         = 0.8,
                       model,      # c("rr_diri_2", "or_diri_2")
                       num.run     = 10,
                       n.adapt     = 2000,
                       n.update    = 5000,
                       n.iter      = 10000,
                       seed = 1){
  library(rjags)

  stopifnot(length(p0) == J)
  stopifnot(length(H0) == J)

  n.dist = 1 + n.exch
  p.dist = rep(1/n.dist, n.dist)
  if (is.null(dist.prior) == TRUE) {dist.prior = rep(1, n.dist)}
  
  tp.index = which(p0 >= H1)
  tn.index = which(p0 <= H0)
  tc.index = which(p0<H1 & p0>H0)
  TP       = NULL
  FP       = NULL
  TN       = NULL
  FN       = NULL
  perfect  = NULL
  FWER     = NULL
  sum.FWPow= 0
  decision = matrix(NA, num.run, J)
  cl.prob  = matrix(NA, num.run, n.dist)
  post.p.dist = matrix(NA, num.run, n.dist)
  post.mean   = NULL
  post.var    = NULL
  
  mat.greater = matrix(NA, num.run, J)

  common = list(J = J, 
                n = n,
                n.dist = n.dist, 
                nex.mean = nex.mean,
                nex.prec = nex.prec,
                ex.mu.mean  = ex.mu.mean,  
                ex.mu.prec  = ex.mu.prec,   
                ex.tau.prec = ex.tau.prec)
  
  for (sim in 1:num.run) {
    set.seed(sim*seed)
    r = rbinom(J, n, p0)

    if (model == "rr_diri_2")    {
      jags.file = "jags_exnex_rr_dirichlet_2.txt"
      sim.data  = append(common, list(r = r, dist.prior = dist.prior))
    }
    if (model == "or_diri_2")    {
      jags.file = "jags_exnex_or_dirichlet_2.txt"
      sim.data  = append(common, list(r = r, dist.prior = dist.prior, H0 = H0))
    }
    
    jags.fit = jags.model(file = jags.file, 
                          data = sim.data,
                          n.adapt = n.adapt, n.chains = 1)
    update(jags.fit, n.burnin = n.update)
    post.sample = jags.samples(jags.fit, 
                               variable.names = c("pi", "cl", "p.dist"),
                               n.iter = n.iter)
    
    pi.chain = post.sample$pi[,,1]
    cl.chain = post.sample$cl[,,1]
    
    post.mean = rbind(post.mean, apply(pi.chain, 1, mean))
    post.var  = rbind(post.var, apply(pi.chain, 1, var))
    
    # traditional decision ----------------------------------------------------
    for (j in 1:J) {
      p.greater  = sum(pi.chain[j,] >= H0[j]) / n.iter
      
      # 0 for nogo, 1 for go
      if (p.greater > tau) {
        decision[sim, j] = 1
      } else {
        decision[sim, j] = 0
      }
      mat.greater[sim, j]  = p.greater
    }
    
    # OC -------------------------------------------------------------------------
    TP[sim] = sum(decision[sim, tp.index] == 1)
    FP[sim] = sum(decision[sim, tn.index] == 1)
    TN[sim] = sum(decision[sim, tn.index] == 0)
    FN[sim] = sum(decision[sim, tp.index] == 0)
    
    perfect[sim] <- (TP[sim] + TN[sim]) == J
    FWER[sim]    <- FP[sim] > 0
    if (length(tp.index)!=0){
      if (TP[sim] == length(tp.index)) {
        sum.FWPow = sum.FWPow + 1
      }
    } 
  }
  
  decision.dist = sapply(1:J, function(j){
    sapply(0:1, function(x){sum(decision[,j]==x) / num.run})
  })
  decision.dist = t(decision.dist)
  colnames(decision.dist) = c("nogo", "go")
  temp = paste0("b", 1:J)
  temp[tp.index] = paste0(temp[tp.index], "(go)")
  temp[tn.index] = paste0(temp[tn.index], "(nogo)")
  temp[tc.index] = paste0(temp[tc.index], "(cons)")
  row.names(decision.dist) = temp

  mean.decision    = colMeans(decision)
  mean.greater  = colMeans(mat.greater)
  
  mean.post.mean = round(apply(post.mean, 2, mean), 4)
  mean.post.var  = round(apply(post.var, 2, mean), 4)
  var.post.mean  = round(apply(post.mean, 2, var),  4)
  mse.post.mean  = round(sapply(1:J, function(j){
    sum((post.mean[,j] - p0[j])^2) / num.run
  }), 4)
  
  out.all = list(cl.chain = cl.chain,  # record of the last simulation
                 decision = decision,
                 mean.greater = mean.greater,
                 rr.post.mean = post.mean,
                 rr.post.var  = post.var)
  out.OC = list(TP = mean(TP),
                FP = mean(FP),
                TN = mean(TN),
                FN = mean(FN),
                perfect = mean(perfect),
                FWER    = mean(FWER),
                FWPR    = sum.FWPow/num.run)
  out.summary = list(mean.post.mean  = mean.post.mean,
                     mean.post.var   = mean.post.var,
                     var.post.mean   = var.post.mean,
                     mse.post.mean   = mse.post.mean,
                     mean.decision   = round(mean.decision,4),
                     mean.greater = round(mean.greater,4),
                     decision.dist   = decision.dist)
  out.setup = list(p0 = p0,
                   model = model,
                   prior = within(sim.data, rm(r)))
  
  out = list(setup  = out.setup,
             mcmc   = out.all,
             summary= out.summary,
             OC     = out.OC,
             readme = "0 for nogo, 1 for go")
  return(out)
}




#########################################################################
p0.from.or <- function(H0, odds_ratio) {    # could be vectors
  p0 = H0*exp(odds_ratio) / (1 - H0 + H0*exp(odds_ratio))
  return(p0)
}



#########################################################################
or.from.p0 <- function(H0, p0) {
  theta = log( (p0/(1-p0)) / (H0/(1-H0)))
  return(theta)
}



