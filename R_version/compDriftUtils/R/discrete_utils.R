#' Simulate two paralogs evolving by compensatory drift.
#'
#' @param trace will plot out the simulations sampled every sample_rate mutations

run_stochastic_sim <- function(p1_start = 5, p2_start = 5, theta = 10, threshold = theta, sim_time = 10^4,
                             mu = 10^-2, N = 10^2, sigma_m = 1, omega = 2,
                             trace = FALSE, sample_rate = 10000, num_sample = 1 ){
  p1_p2 = c()
  for(i in 1:num_sample){
    time <- 0
    p1_0 <- p1_start
    p2_0 <- p2_start
    pp <- 0.5
    if(trace) plot(NULL, xlim = c(0, sim_time), ylim = c(0, 1.1 * theta))
    num_mutations <- 0
    time_trace <- vector("numeric", 10000)
    p1_trace <- vector("numeric", 10000)
    p2_trace <- vector("numeric", 10000)
    index <- 1
    while(time < sim_time){

      ## time of mutation if one or neither is absorbed (abs(p1 - p2) > theta)
      time <- time + rexp(1, rate = 4 * N * mu) #pop. mut. rate for 2 diploid loci; 4N

      if(p1_0 >= threshold | p1_0 <= (theta - threshold)){

        if(p1_0 >= threshold){
          p1_0 <- threshold
          p2_0 <- theta - threshold
        }else{
          p1_0 <- theta - threshold
          p2_0 <- threshold
        }
        break

      }

      if(p2_0 >= threshold | p2_0 <= (theta - threshold)){

        if(p2_0 >= threshold){
          p1_0 <- theta - threshold
          p2_0 <- threshold
        }else{
          p1_0 <- threshold
          p2_0 <- theta - threshold
        }
        break

      }


      if(time > sim_time) break

      p1mut <- p1_0
      p2mut <- p2_0

      if(runif(1) < pp) p1mut <- rnorm(1, p1_0, sigma_m) else p2mut <- rnorm(1, p2_0, sigma_m)
      num_mutations <- num_mutations + 1

      if(runif(1) < pfix(p1_0, p2_0, p1mut, p2mut, theta, omega, N)){
        p1_0 <- p1mut
        p2_0 <- p2mut
      }

      ### store p1 and p2 every X generations
      if((num_mutations %% sample_rate) == 0 & trace){
        lines(cbind(c(time, time_trace[index-1]),
                    c(p1_0, p1_trace[index-1]) ), col = "red")
        lines(cbind(c(time, time_trace[index-1]),
                    c(p2_0, p2_trace[index-1]) ), col = "blue")

        Sys.sleep(0.02)
        time_trace[index] <- time
        p1_trace[index] <- p1_0
        p2_trace[index] <- p2_0
        index <- index + 1
      }
    }#end simulation
    p1_p2 <- rbind(p1_p2, c(p1_0, p2_0))
  }

  colnames(p1_p2) <- c("p1", "p2")
  return(p1_p2)
}

compDrift_var = function(tt){
  return(tt * 1.543 * (mu * omega^3)/(sigma_m * N^(3/2)))
}

fitness = function(xp1, xp2, theta, omega){
  return(exp(-(xp1 + xp2 - theta)^2/(2 * omega^2)))
}

pfix = function(xp1, xp2, xp1_m, xp2_m, theta, omega, N){

  current_fitness = fitness(xp1, xp2, theta, omega)
  mut_fitness = fitness(xp1_m, xp2_m, theta, omega)
  s = mut_fitness/current_fitness - 1
  if(is.nan(s)) s = -1
  fix_prob = (1 - exp(-2 * s))/(1 - exp(-4 * N * s))
  if(s == 0) fix_prob = 1/(2 * N)
  return(fix_prob)

}

