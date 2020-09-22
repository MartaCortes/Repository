library(stats)


insurance <- function(c0, n0, a, tl, lambda, mu, nu, m) {
  capitals <- numeric(m)
  
  for (i in 1:m) {
    t <- 0; n <- n0; capitals[i] <- c0; broken <-  FALSE
    
    while (t < tl && !broken) {
      t_event <- rexp(1, nu + n*mu + n*lambda)
      t <- t + t_event
      
      capitals[i] <- capitals[i] + a*t_event*n 
      
      if (t < tl) {  
        prob <- runif(1)
        if (prob <= nu / (nu + n*mu + n*lambda)) {
          n <- n + 1
        }
        else if (prob <= (n*mu + nu) / (nu + n*mu + n*lambda)) {
          n <- n - 1
        }
        else  {
          capitals[i] <- capitals[i] - rlnorm(1, 5, 1.5) 
          if (capitals[i] < 0) broken = TRUE
        } 
      } else capitals[i] <- capitals[i] + a*(tl-t)*n 
    }
  }
  
  positive_paths <- capitals[capitals>0]
  prop <- length(positive_paths) / length(capitals)
  
  return(list(capital=mean(positive_paths), sd=sd(positive_paths), survival_rate=prop))
}


insurance(c0=10000, n0=50, a=1000, t=20, lambda=2, mu=2, nu = 30, m=10000)
