#' @rdname OU
#' @title Distribution of an Ornstein-Uhlenbeck Process at Time \eqn{t}, Given 
#'   Initial State at Time \eqn{0}
#'   
#' @description An Ornstein-Uhlenbeck (OU) process represents a continuous time 
#'   Markov chain parameterized by an initial state \eqn{x_0}, selection 
#'   strength \eqn{\alpha>0}, long-term mean \eqn{\theta}, and time-unit 
#'   variance \eqn{\sigma^2}. Given \eqn{x_0}, at time \eqn{t}, the state of the
#'   process is characterized by a normal distribution with mean \eqn{x_0 
#'   exp(-\alpha t) + \theta (1 - exp(-\alpha t))} and variance \eqn{\sigma^2 
#'   (1-exp(-2 \alpha t)) / (2 \alpha)}. In the limit \eqn{\alpha -> 0}, the OU 
#'   process converges to a Brownian motion process with initial state \eqn{x_0}
#'   and time-unit variance \eqn{\sigma^2} (at time \eqn{t}, this process is 
#'   characterized by a normal distribution with mean \eqn{x_0} and variance 
#'   \eqn{t \sigma^2}.
#'   
#' @param n Integer, the number of values to sample.
#' @param z Numeric value or vector of size n.
#' @param z0 Numeric value or vector of size n, initial value(s) to condition 
#'   on.
#' @param t Numeric value or vector of size n, denoting the time-step.
#' @param alpha,theta,sigma Numeric values or n-vectors, parameters of the OU 
#'   process; alpha and sigma must be non-negative. A zero alpha is interpreted 
#'   as the Brownian motion process in the limit alpha -> 0.
#' @param log Logical indicating whether the returned density should is on the logarithmic scale.
#'  
#' @details Similar to dnorm and rnorm, the functions described in this
#'   help-page support single values as well as vectors for the parameters z,
#'   z0, t, alpha, theta and sigma.
#' @name OU
NULL

#' @describeIn OU probability density
#' @return dOU returns the conditional probability density(ies) of the elements 
#'   in z, given the initial state(s) z0, time-step(s) t and OU-parameters by
#'   alpha, theta and sigma.
#' @export
dOU <- function(z, z0, t, alpha, theta, sigma, log = TRUE) {
  ett <- exp(-alpha * t)
  sd <- sigma * sqrt(t)
  a <- alpha > 0
  sd[a] <- sigma[a] * sqrt((1 - ett[a]^2) / (2 * alpha[a]))
  mean <- z0 * ett + theta * (1 - ett)
  nan <- is.infinite(t) & alpha == 0
  mean[nan] <- z0[nan]
  
  dnorm(z, mean = mean, sd = sd, log = log)
}

#' @describeIn OU random generator
#'   
#' @return rOU returns a numeric vector of length n, a random sample from the
#'   conditional distribution(s) of one or n OU process(es) given initial
#'   value(s) and time-step(s).
#' @examples 
#' z0 <- 8
#' t <- 10
#' n <- 100000
#' sample <- rOU(n, z0, t, 2, 3, 1)
#' dens <- dOU(sample, z0, t, 2, 3, 1)
#' var(sample)  # around 1/4
#' varOU(t, 2, 1) 
#' 
#' @export
rOU <- function(n, z0, t, alpha, theta, sigma) {
  ett <- exp(-alpha * t)
  sd <- sigma * sqrt(t)
  a <- alpha > 0
  
  sd[a] <- sigma[a] * sqrt((1 - ett[a]^2) / (2 * alpha[a]))
  mean <- z0 * ett + theta * (1 - ett)
  nan <- is.infinite(t) & alpha == 0
  mean[nan] <- z0[nan]
  
  rnorm(n, mean = mean, sd = sd)
}

#' @describeIn OU mean value
#'   
#' @return meanOU returns the expected value of the OU-process at time t.
#' @export
meanOU <- function(z0, t, alpha, theta) {
  ett <- exp(-alpha * t)
  mean <- z0 * ett + theta * (1 - ett)
  nan <- is.infinite(t) & alpha==0
  mean[nan] <- z0[nan]
  mean
}


#' @describeIn OU variance
#'   
#' @return varOU returns the expected variance of the OU-process at time t.
#' @export
varOU <- function(t, alpha, sigma) {
  ett <- exp(-alpha * t)
  var <- sigma^2 * t
  a <- alpha > 0
  var[a] <- sigma[a]^2 * (1 - ett[a]^2) / (2 * alpha[a])
  var
}

#' @describeIn OU standard deviation
#'   
#' @return sdOU returns the standard deviation of the OU-process at time t.
#' @export
sdOU <- function(t, alpha, sigma) {
  ett <- exp(-alpha * t)
  sd <- sigma * sqrt(t)
  a <- alpha > 0
  sd[a] <- sigma[a] * sqrt((1 - ett[a]^2) / (2 * alpha[a]))
  sd
}


#' Generation of a random trajectory of an OU process starting from a given
#' initial state (only for test purpose)
#' @inheritParams rTrajectoryOU
#' 
#' @details Used for test purpose only. This is an internal function and is
#'   appropriate for small time-steps only.
rTrajectoryOUDef <- function(z0, t, alpha, theta, sigma, steps = 1) {
  if(length(t) != steps) {
    t <- rep(t, length.out = steps)
  }
  wien <- sigma * rnorm(steps, 0, sqrt(t))
  z <- numeric(steps + 1)
  z[1] <- z0
  for(i in 1:steps)
    z[i+1] <- z[i] + alpha * (theta - z[i]) * t[i] + wien[i]
  z[-1]
}

#' Generation of a random trajectory of an OU process starting from a given
#' initial state
#' 
#' @description Generates a trajectory xt given initial state z0 according to an
#'   Ornstein-Uhlenbeck process.
#'   
#' @param z0 Numeric value, initial state.
#' @param t Numeric value or vector of size steps, denoting the time-step(s).
#' @param alpha,theta,sigma Numeric values, parameters of the OU process.

#' @param steps Integer, number of steps.
#' @return A numeric vector of length steps containing the generated values
#'   at times 0+t, 0+2t, ..., 0+steps*t.
#'   
#' @examples 
#' z0 <- 8
#' nSteps <- 100
#' t <- 0.01
#' trajectory <- rTrajectoryOU(z0, t, 2, 2, 1, steps = nSteps)
#' plot(trajectory, type = 'l')
#' 
#' @export
rTrajectoryOU <- function(z0, t, alpha, theta, sigma, steps = 1) {
  if(length(t) != steps) {
    t <- rep(t, length.out = steps)
  }
  z <- numeric(steps + 1)
  z[1] <- z0
  ett <- exp(-alpha * t)
  sd <- sigma * sqrt(t)
  a <- alpha > 0
  sd[a] <- sigma[a] * sqrt((1 - ett[a]^2) / (2 * alpha[a]))
  deltas <- rnorm(steps, mean = theta * (1 - ett), sd = sd)
  
  for(i in 1:steps) {
    z[i+1] <- z[i] * ett[i] + deltas[i]
  }
  
  z[-1]
}
