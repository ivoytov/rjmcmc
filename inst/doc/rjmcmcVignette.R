## ------------------------------------------------------------------------
x3 = function(x){
 return(x^3)
}
y = rjmcmc::adiff(x3, c(5,6))
attr(y, "gradient")

## ----echo=FALSE----------------------------------------------------------
load("results")

## ----eval=FALSE----------------------------------------------------------
#  g1 = function(psi){ c(psi[1], log(psi[2]/mu)) }
#  ginv1 = function(theta){ c(theta[1], mu*exp(theta[2])) }

## ----eval=FALSE----------------------------------------------------------
#  g2 = function(psi){ psi }
#  ginv2 = function(theta){ theta }

## ----eval=FALSE----------------------------------------------------------
#  L1 = function(theta){ sum(dpois(y, theta[1], log=T)) }
#  L2 = function(theta){ sum(dnbinom(y, 1/theta[2], mu=theta[1], log=T)) }

## ----eval=FALSE----------------------------------------------------------
#  p.prior1 = function(theta){dgamma(theta[1], lamprior[1], lamprior[2], log=T)
#                             + dnorm(theta[2], 0, sigma, log=T)}
#  

## ----eval=FALSE----------------------------------------------------------
#  p.prior2 = function(theta){dgamma(theta[1], lamprior[1], lamprior[2], log=T)+
#                             dgamma(theta[2], kapprior[1], kapprior[2], log=T)}

## ----eval=FALSE----------------------------------------------------------
#  n = length(y)
#  mu=0.015; sigma=1.5
#  lamprior = c(25,10); kapprior = c(1,10)  # hyperparameters for lambda & kappa
#  
#  goals_post = rjmcmcpost(post.draw = list(draw1, draw2), g = list(g1, g2),
#                 ginv = list(ginv1, ginv2), likelihood = list(L1, L2),
#                 param.prior = list(p.prior1, p.prior2),
#                 model.prior = c(0.5, 0.5), chainlength = 10000)

## ------------------------------------------------------------------------
goals_post

## ----eval=FALSE----------------------------------------------------------
#  dgomp = function(t, A, b, c){ A*exp(-b*exp(-c*t)) }
#  dbert = function(t, L, t0, k){ L*(1-exp(-k*(t+t0))) }

## ----eval=FALSE----------------------------------------------------------
#  L1 = function(theta){sum(dnorm(y, dgomp(t, theta[1], theta[2], theta[3]),
#                                 1/sqrt(theta[4]), log=TRUE))}
#  L2 = function(theta){sum(dnorm(y, dbert(t, theta[1], theta[2], theta[3]),
#                                 1/sqrt(theta[4]), log=TRUE))}

## ----eval=FALSE----------------------------------------------------------
#  g2 = function(psi){ psi }
#  ginv2 = function(theta){ theta }

## ----eval=FALSE----------------------------------------------------------
#  g1 = function(psi){
#    temp = exp(-psi[2]*psi[3])
#    c(psi[1],
#      -log(1-temp),
#      -log((log(1-temp*exp(-psi[3]*tstar))) / (log(1-temp)))/tstar,
#      psi[4])
#  }
#  ginv1 = function(theta){
#    temp = -log((exp(-theta[2]*exp(-theta[3]*tstar))-1)
#                / (exp(-theta[2])-1))/tstar
#    c(theta[1],
#      -log(1-exp(-theta[2]))/temp,
#      temp,
#      theta[4])
#  }

## ----eval=FALSE----------------------------------------------------------
#  p.prior1 = function(theta){
#    sum(dnorm(theta[1:3], 0, 1/sqrt(c(1e-6, 0.05, 1)), log=T)) +
#      dgamma(theta[4], 0.01, 0.01, log=T)
#  }
#  

## ----eval=FALSE----------------------------------------------------------
#  library("FSAdata")
#  data("Croaker2")
#  CroakerM = Croaker2[which(Croaker2$sex=="M"),]
#  y = CroakerM$tl; t = CroakerM$age
#  n = length(y)
#  tstar = 6      # chosen for algorithmic efficiency
#  
#  growth_post=rjmcmcpost(post.draw = list(draw1,draw2), g = list(g1,g2),
#                 ginv = list(ginv1,ginv2), likelihood = list(L1,L2),
#                 param.prior = list(p.prior1,p.prior1),
#                 model.prior = c(0.5,0.5), chainlength = 1e6)

## ----echo=FALSE----------------------------------------------------------
growth_post

## ----eval=FALSE----------------------------------------------------------
#  library("R2jags")
#  
#  inits = function(){list("lambda" = rgamma(1, 1, 0.1),
#                          "kappa" = rgamma(1, 1, 0.1))}
#  params = c("lambda", "kappa")
#  
#  jagsfit1 = jags(data = c('y', 'n', 'lamprior', 'sigma'), inits, params,
#                  n.iter=10000, model.file = "goalsPois.txt")
#  
#  jagsfit2 = jags(data = c('y', 'n', 'lamprior', 'kapprior'), inits, params,
#                  n.iter=10000, model.file = "goalsNB.txt")

## ----eval=FALSE----------------------------------------------------------
#  # Manually
#  fit1 = as.mcmc(jagsfit1); C1 = as.matrix(fit1)
#  draw1 = function(){rev(C1[sample(dim(C1)[1], 1, replace=T),
#                            -which(colnames(C1) == "deviance")])}
#  
#  # Using getsampler function
#  getsampler(jagsfit2, "draw2", order=c(3,2))   # alphabetically, lambda is 3rd parameter and kappa is 2nd

## ----eval=FALSE----------------------------------------------------------
#  ## Gompertz model
#  inits = function(){list(A = abs(rnorm(1, 350, 200)), b = abs(rnorm(1, 2, 3)),
#                          c = abs(rnorm(1, 1, 2)), tau = rgamma(1, 0.1, 0.1))}
#  params = c("A", "b", "c", "tau")
#  jagsfit1 = jags(data = c('y', 't', 'n'), inits, params, n.iter=1e5,
#                  n.thin=20, model.file = "fishGomp.txt")
#  
#  ## von Bertalanffy model
#  inits = function(){list(L = abs(rnorm(1, 350, 200)), t0 = abs(rnorm(1, 2, 3)),
#                          k = abs(rnorm(1, 1, 2)), tau = rgamma(1, 0.1, 0.1))}
#  params = c("L", "t0", "k", "tau")
#  jagsfit2 = jags(data = c('y', 't', 'n'), inits, params, n.iter=1e5,
#                  n.thin=20, model.file = "fishBert.txt")
#  
#  ## Define samplers
#  getsampler(jagsfit1, "draw1")
#  getsampler(jagsfit2, "draw2", c(3,4,2,5))

