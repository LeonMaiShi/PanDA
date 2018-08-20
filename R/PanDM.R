## Script for PanDM Package V1.0.1 ##
## Mar-10-2017 by Leon SHI ##
PanDM <- function(dat, K, epsilon = 1e-04, max.iter = 5000, save.initial = TRUE, para.cores = 0){

  # initialize
  xrow <- nrow(dat); xcol <- ncol(dat)
  loglike0 = loglike1 = list()
  p <- rep(1, K)/K
  q <- matrix(runif(K * xcol), K, xcol)
  mu = sigma = matrix(0, 2, xcol) # Note: here the first row is DM
  mu[1,] <- colMeans(dat, na.rm = TRUE)
  mu[2,] <- mu[1,] + runif(1, 1, 2)
  sigma[1,] <- colSds(dat, na.rm = TRUE)
  sigma[2,] <- sigma[1,] + runif(1)
  if(save.initial)save(p, q, mu, sigma, file = 'Initial_Parameters.rda')

  # function for computing log-likelihood
  myf <- function(a, b, c){
    d <- dnorm(a, mean = b, sd = c, log = TRUE)
    result <- as.vector(d)
    result[is.na(result)] <- 0
    result
  }

  # function for computing posterior cluster membership
  cal_posmem <- function(i, j){
    templike <- matrix(0, xrow, 2) # Note: here the first row is non-DM
    templike[,1] <- log(1 - q[j,i]) + loglike0[[i]]
    templike[,2] <- log(q[j,i]) + loglike1[[i]]
    tempmax <- rowMaxs(templike)
    templike <- exp(templike - tempmax)
    tempsum <- rowSums(templike)
    a <- tempmax + log(tempsum)
    b <- templike / tempsum
    return(list(clustlike = a, condlike =b))
  }

  # EM algorithm to get MLE of p and q
  condlike = list()
  loglike.old <- -1e10

  for(i.iter in 1:max.iter){
	if(i.iter %% 50 == 0 & i.iter < max.iter){
		print(paste0("The first ", i.iter, " iterations have been finished for K=", K, "!"))
	} else if(i.iter == max.iter)print(paste0("The total ", i.iter, " iterations have been finished for K=", K, "!"))

    err <- epsilon + 1

    # compute log-likelihood
    for(i in 1:xcol){
      loglike1[[i]] <- myf(dat[,i], mu[1,i], sigma[1,i])
      loglike0[[i]] <- myf(dat[,i], mu[2,i], sigma[2,i])
    }

    # parallelly compute posterior cluster membership
    clustlike <- matrix(0, xrow, K)

    if(para.cores == 0){
      for(j in 1:K){
        temp_all <- lapply(1:xcol, cal_posmem, j)
        clustlike[,j] <- Reduce('+', lapply(temp_all, function(x)return(x$clustlike))) + log(p[j])
        condlike[[j]] <- do.call('cbind', lapply(temp_all, function(x)return(x$condlike[,2])))
      }
    } else {
      for(j in 1:K){
        temp_all <- mclapply(1:xcol, cal_posmem, j, mc.cores = para.cores)
        clustlike[,j] <- Reduce('+', lapply(temp_all, function(x)return(x$clustlike))) + log(p[j])
        condlike[[j]] <- do.call('cbind', lapply(temp_all, function(x)return(x$condlike[,2])))
      }
    }

    tempmax <- rowMaxs(clustlike)
    clustlike <- exp(clustlike - tempmax)
    tempsum <- rowSums(clustlike)

    # update cluster occurrence rate
    clustlike <- clustlike / tempsum
    clustpsum <- colSums(clustlike)
    p.new <- clustpsum / xrow

    # update cluster patterns
    q.new = matrix(0, K, xcol)
    p.post_k <- list()
    for(j in 1:K){
      p.post_k[[j]] <- clustlike[,j] * condlike[[j]]
      q.new[j,] <- colSums(p.post_k[[j]]) / clustpsum[j]
    }

    # compute posterior p
    p.post <- Reduce('+', p.post_k)

    # update mean and sigma
    mu.new = sigma.new = mu
    mu.new[1,] <- colSums(dat * p.post, na.rm = TRUE) / colSums((p.post * dat / dat), na.rm = TRUE)
    mu.new[2,] <- colSums(dat * (1 - p.post), na.rm = TRUE) / colSums((1 - p.post) * dat / dat, na.rm = TRUE)
    tempmu1 <- matrix(rep(mu.new[1,], xrow), xrow, byrow = TRUE)
    tempmu2 <- matrix(rep(mu.new[2,], xrow), xrow, byrow = TRUE)
    sigma.new[1,] <- sqrt(colSums((dat - tempmu1)^2 * p.post, na.rm = TRUE) / colSums((p.post * dat / dat), na.rm = TRUE))
    sigma.new[2,] <- sqrt(colSums((dat - tempmu2)^2 * (1 - p.post), na.rm = TRUE) / colSums((1 - p.post) * dat / dat, na.rm = TRUE))

    if(all(!is.na(q.new), !is.na(mu.new + sigma.new))){

      # evaluate convergence
      err.p <- max(abs(p.new - p))
      err.q <- max(abs(q.new - q))
      err.mu <- max(abs(mu.new - mu))
      err.sigma <- max(abs(sigma.new - sigma))
      err <- max(err.p, err.q, err.mu, err.sigma)

      # evaluate whether the log.likelihood increases
      loglike.new <- sum(tempmax + log(tempsum))
      delta_loglike <- loglike.new - loglike.old

      # update parameters
      p <- p.new
      q <- q.new
      mu <- mu.new
      sigma <- sigma.new
      loglike.old <- loglike.new

      if(err < epsilon){
        print(paste0("PanDM ended as the maximum absolute error is smaller than ", epsilon, "!"))
        break;
      }
    }
    else{
      print("PanDM stopped due to the appearance of NaN!")
      break;
    }
  }

  # compute posterior p (i.e. Pr(H_{gc} = 1))
  for(i in 1:xcol) {
    loglike1[[i]] <- myf(dat[,i], mu[1,i], sigma[1,i])
    loglike0[[i]] <- myf(dat[,i], mu[2,i], sigma[2,i])
  }

  if(para.cores == 0){
    for(j in 1:K){
      temp_all <- lapply(1:xcol, cal_posmem, j)
      clustlike[,j] <- Reduce('+', lapply(temp_all, function(x)return(x$clustlike))) + log(p[j])
      condlike[[j]] <- do.call('cbind', lapply(temp_all, function(x)return(x$condlike[,2])))
    }
  } else {
    for(j in 1:K){
      temp_all <- mclapply(1:xcol, cal_posmem, j, mc.cores = para.cores)
      clustlike[,j] <- Reduce('+', lapply(temp_all, function(x)return(x$clustlike))) + log(p[j])
      condlike[[j]] <- do.call('cbind', lapply(temp_all, function(x)return(x$condlike[,2])))
    }
  }

  tempmax <- rowMaxs(clustlike)
  clustlike <- exp(clustlike - tempmax)
  tempsum <- rowSums(clustlike)
  clustlike <- clustlike/tempsum

  for(j in 1:K)p.post_k[[j]] <- clustlike[,j] * condlike[[j]]

  p.post <- Reduce('+', p.post_k)

  cluster <- apply(clustlike, 1, which.max)

  # calculate BIC
  bic <- -2 * loglike.old + (K - 1 + K * xcol + xcol * 2 * 2) * log(xrow)

  # results
  result <- list(cluster = cluster, p.post = p.post, p = p, q = q, loglike = loglike.old, bic = bic, mu = mu, sigma = sigma)
}

PanDM_BIC <- function(resultlist, BICplot = TRUE){
  myK <- sapply(resultlist, function(x)return(length(x$p)))
  bictable <- cbind(myK, sapply(resultlist, function(x)return(x$bic)))
  colnames(bictable) <- c('K','BIC')
  print(bictable)

  if(BICplot){
    K_i <- which.min(bictable[,2]);
    K <- bictable[K_i,1] # record the K with minimal BIC

    plot(bictable[,2]~bictable[,1], pch = 16, xlab = "Number of Patterns (K)", cex = 1.5,
         ylab = "Bayesian Information Criterion", main = "BIC for PanDM model fitting",
         xaxt = 'none', cex.main = 1, cex.lab = 1, cex.axis = 1)
    lines(bictable[,2]~bictable[,1], lwd=2)
    axis(1, at = myK, cex.axis = 1)
    points(x = K, y = bictable[K_i,2], pch = 11, col = 'red', bg = 'red', lwd = 3)
  }
}
