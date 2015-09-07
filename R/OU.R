#' Random generation of values given an initial value following an OU process.
#' 
#' @description
#' Generates xt given x0 according to an Ornstein-Uhlenbeck process with 
#' parameters alpha, theta and sigma. Uses the formula based on a scaled 
#' time transformed Wiener process. See
#' http://en.wikipedia.org/wiki/Ornstein–Uhlenbeck_process#Alternative_representation_for_nonstationary_processes
#' 
#' @details Unlike the function OU.def this one works well with big time-steps 
#' as well.
#' 
#' @param x0 initial value
#' @param dt time step. Doesn't need to be infinitesimally small
#' @param alpha the strength of the selective constraint
#' @param theta long term mean value of the OU process
#' @param sigma the unit-time standard deviation of the random component in the OU
#'     process.
#' @param steps number of steps of the OU process
#
#' @return
#'   A numeric vector of length steps containing the values of the OU process 
#'     at times 0+dt, 0+2dt, ..., 0+steps*dt.
#' @export
OU <- function(x0, dt, alpha, theta, sigma, steps=1) {
  x <- numeric(steps+1)
  x[1] <- x0
  if(alpha == 0) {
    deltas <- sigma*rnorm(steps, 0, sqrt(dt))
  } else {
    # alpha >0
    deltas <- theta*(1-exp(-alpha*dt)) + 
      sigma/sqrt(2*alpha)*exp(-alpha*dt)*rnorm(steps, 0, sqrt(exp(2*alpha*dt)-1))
  }
  
  for(t in 1:steps)
    x[t+1] <- x[t]*exp(-alpha*dt)+deltas[t]
  x[-1]
}

#' Random generation of values given an initial value following an OU process. 
#' @details
#' Used for test purpose only. This is an internal function and is appropriate 
#' for small time-steps only.
#' 
#' @description
#' Generates xt given x0 according to an Ornstein-Uhlenbeck process with 
#' parameters theta, alpha and sigma. Works fine for small enough dt. 
OU.def <- function(x0, dt, alpha, theta, sigma, steps=1) {
  wien <- sigma*rnorm(steps, 0, sqrt(dt))
  x <- numeric(steps+1)
  x[1] <- x0
  for(t in 1:steps)
    x[t+1] <- x[t] + alpha*(theta-x[t])*dt + wien[t]
  x[-1]
}

#' non-incremental version of the OU process based on the scaled time-
#' transformed Wiener process formula. For test purposes only.
OU.noninc <- function(x0, alpha, theta, sigma, from, to, steps=2) {
  t<-seq(from, to, length=steps)
  x0*exp(-alpha*t) + theta*(1-exp(-alpha*t)) + 
    rnorm(n=steps, mean=0, sd=sqrt(sigma^2/(2*alpha)*(1-exp(-2*alpha*t))))
}

#' Conditional probability density of a value following an OU process, given 
#' initial value and time-step
#' @export
dCondOU <- function(xt, x0, t, alpha, theta, sigma, log=T) {
  ett <- exp(-alpha*t)
  if(alpha == 0) {
    dnorm(xt, mean=x0, sd=sigma*sqrt(t), log=log)
  } else {
    # alpha > 0
    dnorm(xt, mean=x0*ett+theta*(1-ett), sd=sqrt(sigma^2/(2*alpha)*(1-ett^2)), 
          log=log)
  }
}

#' Probability density of the stationary distribution of an OU process.
#' @export
dStOU <- function(x, alpha, theta, sigma, log=T) {
  dnorm(x, theta, sigma/sqrt(2*alpha), log=log)
}

#' Random sampling from the stationary distribution of an OU process.
#' @export
rStOU <- function(n, alpha, theta, sigma) {
  rnorm(n, theta, sigma/sqrt(2*alpha))
}

#' Random sampling from the conditional distribution of an OU process given 
#' initial value and time-step
#' @export
rCondOU <- function(n, x0, t, alpha, theta, sigma) {
  ett <- exp(-alpha*t)
  rnorm(n, mean=x0*ett+theta*(1-ett), sd=sqrt(sigma^2/(2*alpha)*(1-ett^2)))
}

#' Variance of the stationary distribution of an OU-process.
#' @export
varStOU <- function(alpha, sigma) {
  sigma^2/(2*alpha)
}

#' Standard deviation of the stationary distribution of an OU-process
#' @export
sdStOU <- function(alpha, sigma) {
  sigma/sqrt(2*alpha)
}