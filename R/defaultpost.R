
#' Perform Post-Processing Using Default Bijections
#' 
#' Performs Bayesian multimodel inference, estimating Bayes factors and 
#' posterior model probabilities for N candidate models. Unlike 
#' \code{\link{rjmcmcpost}}, this function uses a default bijection scheme based
#' on approximating each posterior by a multivariate normal distribution. The 
#' result is reminiscent of the algorithm of Carlin & Chib (1995) with a 
#' multivariate normal pseudo-prior. Transformation Jacobians are computed using
#' automatic differentiation so do not need to be specified.
#' 
#' @param coda A list of N codas corresponding to the posterior distribution 
#'   under each model. Generally this coda output will be found using MCMC. Note
#'   that each parameter should be real-valued so some parameters may need to be
#'   transformed, using logarithms for example.
#' @param likelihood A list of N functions specifying the log-likelihood 
#'   functions for the data under each model.
#' @param param.prior A list of N functions specifying the prior distributions 
#'   for each model-specific parameter vector.
#' @param model.prior A numeric vector of the prior model probabilities. Note 
#'   that this argument is not required to sum to one as it is automatically 
#'   normalised.
#' @param chainlength How many iterations to run the Markov chain for.
#' @param TM.thin How regularly to calculate transition matrices as the chain 
#'   progresses.
#' @param progress A logical determining whether a progress bar is drawn.
#' @return Returns a list-like \code{rj} object with named elements 
#'   \code{result}, \code{densities}, \code{psidraws},
#'   \code{progress} and \code{meta}. See \code{\link{rjmethods}} for
#'   information on what each of these elements are.
#' 
#' @importFrom stats var  
#' @references Carlin, B. P. and Chib, S. (1995) Bayesian Model Choice via 
#'   Markov Chain Monte Carlo Methods. \emph{Journal of the Royal Statistical 
#'   Society, Series B, 473-484}.
#' @references Barker, R. J. and Link, W. A. (2013) Bayesian multimodel 
#'   inference by RJMCMC: A Gibbs sampling approach. \emph{The American 
#'   Statistician, 67(3), 150-156}.
#'   
#' @seealso \code{\link{adiff}} \code{\link{rjmcmcpost}}
#'   
#' @examples
#' ## Comparing two binomial models -- see Barker & Link (2013) for further details.
#' 
#' y=c(8,16); sumy=sum(y)
#' n=c(20,30); sumn=sum(n)
#' 
#' L1=function(p){if((all(p>=0))&&(all(p<=1))) sum(dbinom(y,n,p,log=TRUE)) else -Inf}
#' L2=function(p){if((p[1]>=0)&&(p[1]<=1)) sum(dbinom(y,n,p[1],log=TRUE)) else -Inf}
#' 
#' p.prior1=function(p){sum(dbeta(p,1,1,log=TRUE))}
#' p.prior2=function(p){dbeta(p[1],1,1,log=TRUE)+dbeta(p[2],17,15,log=TRUE)}
#' 
#' coda1=matrix(rbeta(2000,y+1,n-y+1), 1000, 2, byrow=TRUE)  ## full conditional posterior
#' coda2=matrix(c(rbeta(1000,sumy+1,sumn-sumy+1),rbeta(1000,17,15)), 1000, 2)
#' 
#' out=defaultpost(coda=list(coda1,coda2), likelihood=list(L1,L2), 
#'                 param.prior=list(p.prior1,p.prior2), model.prior=c(1,1), chainlength=1000)
#' 
#' @export
defaultpost=function(coda, likelihood, param.prior, model.prior, chainlength=10000, TM.thin=chainlength/10, progress=TRUE){
  n.models = length(coda)
  nTM = chainlength/TM.thin
  TM = rep(list(matrix(NA, n.models, n.models)), nTM); mvnd = rep(NA, n.models)
  
  n.par = rep(NA, n.models)
  for(j in 1:n.models){
    coda[[j]] = as.matrix(coda[[j]])
    if(any(colnames(coda[[j]])=="deviance")){ 
      coda[[j]] = coda[[j]][, -which(colnames(coda[[j]])=="deviance")]
    }
    n.par[j] = ncol(coda[[j]])
  }
  dim.psi = sum(n.par)
  
  p.bar = psi = rep(NA, dim.psi)
  modlab = c(); covar = list()
  for(j in 1:n.models){
    modlab = c(modlab, rep(j, n.par[j]))
    post = coda[[j]]
    p.bar[which(modlab==j)] = apply(post, 2, mean)
    covar[[j]] = var(post)  # covariance matrix for each model
  }
  u.prior = cbind(modlab, p.bar)
  
  psistore = matrix(NA, chainlength, dim.psi)
  store = rep(list(matrix(NA, chainlength, n.models*3, dimnames=list(NULL, c(paste0("Posterior M", 1:n.models), paste0("Likelihood M", 1:n.models), paste0("Prior M", 1:n.models))))), n.models)
  message('Post-Processing Based on Normal Pseudo-Prior')
  for(j in 1:n.models){
    message('Row ', j, appendLF=FALSE)
    wuse = trunc(getOption("width")-20L)    # Set up progress bar
    if(progress){ pb = utils::txtProgressBar(min=0, max=chainlength, initial=0, char="*", style=3, width=wuse) }
    
    term = matrix(NA,chainlength,n.models)
    is = which(modlab==j)
    
    for(i in 1:chainlength){   
      psi[is] = coda[[j]][sample(dim(coda[[j]])[1], 1),]
      for(m in 1:n.models){
        if(m==j){ next }
        mis = which(modlab==m)
        psi[mis] = mvtnorm::rmvnorm(1, u.prior[mis,2], covar[[m]])
      }
      psistore[i,] = psi
      
      for(k in 1:n.models){
        ind = which(modlab==k)
        like = likelihood[[k]](psi[ind])
        prior = param.prior[[k]](psi[ind])
        for(m in 1:n.models){
          if(m==k){ next }
          mis = which(modlab==m)
          prior = prior + sum(mvtnorm::dmvnorm(psi[mis], u.prior[mis,2], covar[[m]], log=T))  
        }
        term[i,k] = like + prior + log(model.prior[k])
        store[[j]][i, k+n.models*(0:2)] = c(term[i,k], like, prior)
      }
      term[i,] = term[i,] - max(term[i,])
      term[i,] = exp(term[i,])/sum(exp(term[i,]))
      if(any(is.na(term[i,]))){ warning(paste("NAs in chain for model",j)); break }
      if(progress){ utils::setTxtProgressBar(pb, value=i) }
      if(i%%TM.thin == 0){
        TM[[i/TM.thin]][j,]=apply(term[1:i,], 2, mean)
      }
    }
    if(progress){ close(pb) }
  }
  
  prob = BF = matrix(NA,nTM,n.models)
  for(i in 1:nTM){
    ev = eigen(t(TM[[i]]))
    prob.us = ev$vector[,which(abs(ev$values-1) < 1e-8)]
    prob[i,] = prob.us/sum(prob.us)
    BF[i,] = prob[i,]/prob[i,1] * model.prior[1]/model.prior
  }
  return(rj(list(result=list("Transition Matrix" = TM[[nTM]], "Posterior Model Probabilities"=prob[nTM,], 
                             "Bayes Factors" = BF[nTM,], "Second Eigenvalue" = ev$value[2]), 
                 densities = store, psidraws = psistore, progress=list(TM=TM, prb=prob), 
                 meta=list(chainlength=chainlength, TM.thin=TM.thin))))
}
