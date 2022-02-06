pitaun <- function(n_max, theta){
  nrows <- length(theta)  
  mu <- cos(theta)
  # Initialize recurrence pi_0=0, pi_1=1
  # pinm1 contains pi_{n-1} n=1..nn_max+1
  pinm1 <- matrix(0, nrows, n_max + 1)
  pinm1[ , 2] <- 1
  # Get pi_2 to pi_nn_max by recurrence (see SPlaC guide)
  # pi_n is pinm1(:,n+1)
  for (n in seq(2, n_max)) {
    pinm1[ ,n+1] <- (2*n-1)/(n-1) * mu * pinm1[ ,n]- n/(n-1)*pinm1[,n-1]
  }
  pin <- pinm1[,-1L]
  # return tau_n matrix
  nvec <- seq_len(n_max)
  nmat <- matrix(nvec+1, ncol=n_max, nrow=nrows, byrow = TRUE) 
  taun <- outer(mu,nvec) * pin - nmat * pinm1[,nvec]
  
  # return pi_n matrix (except n=0)
  list(pin = pin, taun = taun)
}
