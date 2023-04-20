


###################################################
##### Get posterior samples from EXNEX Method #####
###################################################



jags.exnex.post.samples.3ODM <- function(nk = rep(24, 6),
                                         rk = rep(4, 6),
                                         K = 6,
                                         Nexch = 2,
                                         Nmix = 3,
                                         pMix = c(0.5, 0, 0.5),
                                         dirichprior = c(1, 1, 1),
                                         dirichprior.check = TRUE,
                                         mu.mean = c(-1.735, 0.847),
                                         mu.prec = c(0.146, 0.266),
                                         tau.HN.scale = c(1, 1),
                                         nex.mean = rep(-1.734, 6),
                                         nex.prec = rep(0.128, 6),
                                         n.burnin = 5000,
                                         n.chain = 1,
                                         n.iter = 8000,
                                         LOR = FALSE,
                                         H0 = NULL) {
  if (dirichprior.check == TRUE & LOR == TRUE) {
    jags.data <- list(
      J = K,
      n = nk,
      r = rk,
      n.dist = Nmix,
      nex.mean = nex.mean,
      nex.prec = nex.prec,
      ex.mu.mean  = mu.mean,
      ex.mu.prec  = mu.prec,
      ex.tau.prec = tau.HN.scale,
      H0 = H0,
      dist.prior = dirichprior
    )
    
    jags.fit <-
      jags.model(
        file = "jags_exnex_or_dirichlet_2.txt",
        data = jags.data,
        n.adapt = 5000,
        n.chains = n.chain,
        quiet = T,
        init = list(".RNG.name" = "base::Wichmann-Hill",
                    ".RNG.seed" = 31415)
      )
    update(jags.fit, n.burnin)
  } else if (dirichprior.check == FALSE & LOR == TRUE) {
    jags.data <- list(
      J = K,
      n = nk,
      r = rk,
      n.dist = Nmix,
      nex.mean = nex.mean,
      nex.prec = nex.prec,
      ex.mu.mean  = mu.mean,
      ex.mu.prec  = mu.prec,
      ex.tau.prec = tau.HN.scale,
      H0 = H0,
      p.dist = pMix
    )
    
    jags.fit <-
      jags.model(
        file = "jags_exnex_or_no_dirichlet.txt",
        data = jags.data,
        n.adapt = 5000,
        n.chains = n.chain,
        quiet = T,
        init = list(".RNG.name" = "base::Wichmann-Hill",
                    ".RNG.seed" = 31415)
      )
    update(jags.fit, n.burnin)
  } else if (dirichprior.check == TRUE & LOR == FALSE) {
    jags.data <- list(
      J = K,
      n = nk,
      r = rk,
      n.dist = Nmix,
      nex.mean = nex.mean,
      nex.prec = nex.prec,
      ex.mu.mean  = mu.mean,
      ex.mu.prec  = mu.prec,
      ex.tau.prec = tau.HN.scale,
      dist.prior = dirichprior
    )
    
    jags.fit <-
      jags.model(
        file = "jags_exnex_rr_dirichlet_2.txt",
        data = jags.data,
        n.adapt = 5000,
        n.chains = n.chain,
        quiet = T,
        init = list(".RNG.name" = "base::Wichmann-Hill",
                    ".RNG.seed" = 31415)
      )
    update(jags.fit, n.burnin)
  } else if (dirichprior.check == FALSE & LOR == FALSE) {
    jags.data <- list(
      J = K,
      n = nk,
      r = rk,
      n.dist = Nmix,
      nex.mean = nex.mean,
      nex.prec = nex.prec,
      ex.mu.mean  = mu.mean,
      ex.mu.prec  = mu.prec,
      ex.tau.prec = tau.HN.scale,
      p.dist = pMix
    )
    
    jags.fit <-
      jags.model(
        file = "jags_exnex_rr_no_dirichlet.txt",
        data = jags.data,
        n.adapt = 5000,
        n.chains = n.chain,
        quiet = T,
        init = list(".RNG.name" = "base::Wichmann-Hill",
                    ".RNG.seed" = 31415)
      )
    update(jags.fit, n.burnin)
  }
  
  post.samples <-
    coda.samples(
      jags.fit,
      variable.names = c("pi"),
      n.iter = n.iter / n.chain,
      progress.bar = "none"
    )
  exnex.out <- do.call(rbind, post.samples)
  
  return(exnex.out)
  
}

#######################################
##### Simulation for EXNEX Method #####
#######################################

exnex.3ODM.sim <- function(num.sim = 1000,
                           p0 = rep(0.2, 6),
                           H0 = rep(0.2, 6),
                           H1 = rep(0.4, 6),

                           K = 6,
                           Nexch = 2,
                           Nmix = 3,
                           pMix = c(0.5, 0, 0.5),
                           dirichprior = c(1, 1, 1),
                           dirichprior.check = TRUE,
                           mu.mean = c(-1.735, 0.847),
                           mu.prec = c(0.146, 0.266),
                           tau.HN.scale = c(1, 1),
                           nex.mean = rep(-1.734, 6),
                           nex.prec = rep(0.128, 6),
                           n.burnin = 4000,
                           n.chain = 1,
                           n.iter = 10000,
                           seed = 1,
                           ###3ODM option
                           # Specify if it's LOR or RR model
                           LOR = FALSE,
                           tau.TV      = 0.2,
                           tau.LRV     = 0.8,
                           #interim analysis criteria
                           L = 0,
                           each.plan   = matrix(rep(15, 6), 1, 6, byrow = T),
                           nogo_LRV    = 0.3,
                           nogo_TV     = 0.05,
                           go_LRV      = 0.95,
                           go_TV       = 0.3) {
  tp <- which(p0 >= H1)
  tn <- which(p0 <= H0)
  tc <- which(p0 < H1 & p0 > H0)
  
  TP <- NULL
  FP <- NULL
  TC <- NULL
  FC <- NULL
  TN <- NULL
  FN <- NULL
  perfect <- NULL
  FWER <- NULL
  decision.final <- matrix(NA, num.sim, K)
  decision.interim <- matrix(NA, num.sim, K)
  decision.fulldata <- matrix(NA, num.sim, K)
  trace.enroll = matrix(NA, num.sim, K)
  
  mat.inconsist = matrix(0, num.sim, 6)
  colnames(mat.inconsist) = c("0->1", "0->2", "1->0", "1->2", "2->0", "2->1")
  
  post.mean = NULL
  post.var  = NULL
  sample.size = matrix(NA, num.sim, K)
  
  mat.greater.TV  = matrix(NA, num.sim, K)
  mat.greater.LRV = matrix(NA, num.sim, K)
  #rk_all<-list()
  
  inter.decision.dist = final.decision.dist = mat.inconsit.dist = NA
  
  for (sim in 1:num.sim) {
    set.seed(sim * seed)
    #new variable stop
    stop <- NULL
    #new variable cont
    cont <- 1:K
    
    each.true     <-
      matrix(NA, L + 1, K)    # truly newly enrolled in each enrollment
    each.response <- matrix(NA, L + 1, K)    #
    #If there is more than one interim analysis
    l = 1
    if (L != 0) {
      dec.interim <- NULL
      
      while (l <= L & length(cont) > 0) {
        dec.interim.vec <- rep(NA, K)
        
        # truly newly enrolled people in each enrollment:-------------------------
        each.true[l, ]     <-
          sapply(1:K, function(j) {
            ifelse(j %in% cont, each.plan[l, j], 0)
          })
        each.response[l, ] <-
          sapply(1:K, function(j) {
            rbinom(1, each.true[l, j], p0[j])
          })
        
        nk = colSums(each.true, na.rm = TRUE)
        rk = colSums(each.response, na.rm = TRUE)
        
        post.sample <-
          jags.exnex.post.samples.3ODM(
            nk,
            rk,
            K,
            Nexch,
            Nmix,
            pMix,
            dirichprior,
            dirichprior.check,
            mu.mean,
            mu.prec,
            tau.HN.scale,
            nex.mean,
            nex.prec,
            n.burnin,
            n.chain,
            n.iter,
            LOR,
            H0
          )
        
        post.sample <- t(post.sample)
        
        # 3ODM --------------------------------------------------------------------
        for (j in 1:K) {
          p.greater.TV  = sum(post.sample[j, ] > H1[j]) / n.iter
          p.greater.LRV = sum(post.sample[j, ] > H0[j]) / n.iter
          
          if (j %in% stop) {
            dec.interim.vec[j] = 0
          } else if (p.greater.LRV < nogo_LRV &
                     p.greater.TV  < nogo_TV) {
            stop = c(stop, j)
            cont = setdiff(cont, stop)
            dec.interim.vec[j] = 0
            # "-n" means an indication is stopped at n-th enrollment:
            trace.enroll[sim, j] <- -l
          } else if (p.greater.LRV > go_LRV &
                     p.greater.TV  > go_TV) {
            dec.interim.vec[j] = 2
          } else{
            dec.interim.vec[j] = 1
          }
        }
        dec.interim <- rbind(dec.interim, dec.interim.vec)
        l = l + 1
      }
      
      M = nrow(dec.interim)
      
      decision.interim[sim, ] = sapply(1:K, function(j) {
        ifelse(max(dec.interim[, j]) == 2, 2, dec.interim[M, j])
      })
    }
    
    #final analysis
    each.true[L + 1, ]     <-
      sapply(1:K, function(j) {
        ifelse(j %in% cont, each.plan[L + 1, j], 0)
      })
    each.response[L + 1, ] <-
      sapply(1:K, function(j) {
        rbinom(1, each.true[L + 1, j], p0[j])
      })
    
    nk = colSums(each.true, na.rm = TRUE)
    rk = colSums(each.response, na.rm = TRUE)
    
    post.sample <-
      jags.exnex.post.samples.3ODM(
        nk,
        rk,
        K,
        Nexch,
        Nmix,
        pMix,
        dirichprior,
        dirichprior.check,
        mu.mean,
        mu.prec,
        tau.HN.scale,
        nex.mean,
        nex.prec,
        n.burnin,
        n.chain,
        n.iter,
        LOR,
        H0
      )
    post.sample <- t(post.sample)
    sample.size[sim, ] = nk
    
    for (j in 1:K) {
      p.greater.TV  = sum(post.sample[j, ] > H1[j]) / n.iter
      p.greater.LRV = sum(post.sample[j, ] > H0[j]) / n.iter
      
      if (p.greater.LRV < tau.LRV & p.greater.TV < tau.TV) {
        decision.fulldata[sim, j] = 0
      } else if (p.greater.LRV > tau.LRV & p.greater.TV > tau.TV) {
        decision.fulldata[sim, j] = 2
      } else {
        decision.fulldata[sim, j] = 1
      }
      mat.greater.TV[sim, j]  = round(p.greater.TV, 4)
      mat.greater.LRV[sim, j] = round(p.greater.LRV, 4)
    }
    
    
    # early nogo in interim will be nogo in final decision
    # whatever decision is in the full-data decision
    if (L >= 1) {
      for (j in 1:K) {
        if (decision.interim[sim, j] == 0) {
          decision.final[sim, j] = 0
        } else {
          decision.final[sim, j] = decision.fulldata[sim, j]
        }
      }
    } else
      (
        decision.final = decision.fulldata
      )
    
    
    for (i in 0:2) {
      for (j in setdiff(0:2,i)) {
        n.incons = sum(decision.final[sim,which(decision.interim[sim,] == i)] == j)
        mat.inconsist[sim, glue::glue("{i}->{j}")] = n.incons
      }
    }
    
    post.mean = rbind(post.mean, apply(post.sample, 1, mean))
    post.var  = rbind(post.var, apply(post.sample, 1, var))
    # OC -------------------------------------------------------------------------
    TP[sim] = sum(decision.final[sim, tp] == 2)
    FP[sim] = sum(decision.final[sim, c(tn, tc)] == 2)
    TN[sim] = sum(decision.final[sim, tn] == 0)
    FN[sim] = sum(decision.final[sim, c(tp, tc)] == 0)
    TC[sim] = sum(decision.final[sim, tc] == 1)
    FC[sim] = sum(decision.final[sim, c(tp, tn)] == 1)
    
    perfect[sim] <- (TP[sim] + TN[sim] + TC[sim]) == K
    FWER[sim]    <- FP[sim] > 0
  }
  
  #Summary
  inter.decision.dist = sapply(1:K, function(j) {
    sapply(0:2, function(x) {
      sum(decision.interim[, j] == x) / num.sim
    })
  })
  
  inter.decision.dist = t(inter.decision.dist)
  colnames(inter.decision.dist) = c("early_nogo", "early_cons", "early_go")
  
  final.decision.dist = sapply(1:K, function(j) {
    sapply(0:2, function(x) {
      sum(decision.final[, j] == x) / num.sim
    })
  })
  final.decision.dist = t(final.decision.dist)
  colnames(final.decision.dist) = c("final_nogo", "final_cons", "final_go")
  
  if (L != 0) {
    mat.inconsist = mat.inconsist[, c(1, 2, 5, 6)]
    mat.inconsit.dist = colMeans(mat.inconsist)
    names(mat.inconsit.dist) = c("NG->Consider", "NG->Go", "Go->NG", "Go->Consider")
  }
  
  out <- list(
    inter.decision.dist = inter.decision.dist,
    final.decision.dist = final.decision.dist,
    mat.inconsist.dist     =   mat.inconsit.dist,
    #             post.mean=post.mean,
    #             mat.greater.TV    = mat.greater.TV,
    #             mat.greater.LRV   = mat.greater.LRV,
    sample.size = apply(sample.size, 2, mean)
    #             TP=mean(TP), FP=mean(FP),TN=mean(TN), FN=mean(FN), perfect=mean(perfect), FWER=mean(FWER)
  )
  return(out)
}

#



test <-
  exnex.3ODM.sim(
    L = 3,
    num.sim = 20,
    p0 = rep(0.1, 6),
    each.plan = matrix(rep(20, 24), 4, 6, byrow = T),
    LOR = TRUE
  )