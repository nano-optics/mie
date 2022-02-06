library(Bessel)
library(BH)
library(Rcpp)

sourceCpp(code='
          #include <Rcpp.h>
          // [[Rcpp::depends(BH)]]
          #include <boost/math/special_functions/bessel.hpp>
          using namespace boost::math;
          // [[Rcpp::export]]
          double bj(const int nu, const double x) {
          double y = cyl_bessel_j(nu, x);
          return (y);
          }'
)


ricatti_bessel <- function(rho, n_max){
  
  # rho <- as.complex(rho)
  sq <- sqrt((pi/2)*rho)
  # sq <- 1
  bj <- sq * BesselJ(z = rho, nu = 0.5, nSeq = n_max+1)
  bh <- sq * BesselH(z = rho, nu = 0.5, nSeq = n_max+1, m = 1)
  
  psi <- bj[, -1L]
  xi <- bh[, -1L]
  norho <- outer(1/rho, seq_len(n_max))
  
  dpsi <- bj[, -(n_max+1)] - norho * psi
  dxi <- bh[, -(n_max+1)] - norho * xi
  
  list(psi=psi, xi=xi, dpsi=dpsi, dxi=dxi)
}

rho <- 2*pi/seq(500, 800, by=10) * 50 * 1.5
n_max <- 5
r <- ricatti_bessel(rho, n_max)

cpp <- sapply(rho, bj, nu=0.5)
matplot(rho, cbind(cpp, r$psi[,1]), t="l",lty=1)
