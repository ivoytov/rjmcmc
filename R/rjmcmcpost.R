
#' Perform Reversible-Jump MCMC Post-Processing
#' 
#' Performs Bayesian multimodel inference, estimating Bayes factors and 
#' posterior model probabilities for N candidate models. Using the 'universal 
#' parameter' restriction in Barker & Link (2013), RJMCMC is treated as a Gibbs 
#' sampling problem, where the algorithm alternates between updating the model 
#' and the model specific parameters. Transformation Jacobians are computed 
#' using automatic differentiation so do not need to be specified.
#' 
#' @param post.draw A list of N functions that randomly draw from the posterior 
#'   distribution under each model. Generally these functions sample from the 
#'   coda output of a model fitted using MCMC. Functions that draw from the 
#'   posterior in known form are also allowed.
#' @param g A list of N functions specifying the bijections from the universal 
#'   parameter \code{psi} to each model-specific parameter set.
#' @param ginv A list of N functions specifying the bijections from each 
#'   model-specific parameter set to \code{psi}. These are the inverse 
#'   transformations of \code{g}.
#' @param likelihood A list of N functions specifying the log-likelihood 
#'   functions for the data under each model.
#' @param param.prior A list of N functions specifying the prior distributions 
#'   for each model-specific parameter vector.
#' @param model.prior A numeric vector of the prior model probabilities. Note 
#'   that this argument is not required to sum to one as it is automatically 
#'   normalised.
#' @param chainlength How many iterations to run the Markov chain for.
#' @return Returns a list with three elements. The first element, \code{TM}, 
#'   gives the transition matrix for the Markov chain, summarising how the chain
#'   moved over time. The second element, \code{prb}, gives the posterior model 
#'   probability for each model. The third element, \code{BF}, gives the Bayes 
#'   factor in favour of each model compared to the first model.
#'   
#' @references Barker, R. J. and Link, W. A. (2013) Bayesian multimodel 
#'   inference by RJMCMC: A Gibbs sampling approach. \emph{The American 
#'   Statistician, 67(3), 150-156}.
#'   
#' @seealso \code{\link{adiff}} \code{\link{getsampler}}
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
#' g1=function(psi){p=psi}
#' g2=function(psi){w=n[1]/sum(n); p=c(w*psi[1]+(1-w)*psi[2],psi[2])}
#' ginv1=function(p){p}
#' ginv2=function(p){c(sum(n)/n[1]*p[1]-n[2]/n[1]*p[2],p[2])}
#' 
#' p.prior1=function(p){sum(dbeta(p,1,1,log=TRUE))}
#' p.prior2=function(p){dbeta(p[1],1,1,log=TRUE)+dbeta(p[2],17,15,log=TRUE)}
#' 
#' draw1=function(){rbeta(2,y+1,n-y+1)}
#' draw2=function(){c(rbeta(1,sumy+1,sumn-sumy+1),rbeta(1,17,15))}
#' 
#' out=rjmcmcpost(post.draw=list(draw1,draw2), g=list(g1,g2), ginv=list(ginv1,ginv2), 
#'                likelihood=list(L1,L2), param.prior=list(p.prior1,p.prior2), 
#'                model.prior=c(0.5,0.5), chainlength=10000)
#' 
#' @export
rjmcmcpost=function(post.draw, g, ginv, likelihood, param.prior, model.prior, chainlength=10000){
  n.models = length(post.draw)
  TM = matrix(NA,n.models,n.models)
  
  detJtest = matrix(NA, 100, n.models); samedet=rep(FALSE, n.models)
  
  for(j in 1:n.models){
    message('Row ',j,appendLF=FALSE)    # Set up progress bar
    wuse = trunc(getOption("width")-20L)
    pb = utils::txtProgressBar(min=0,max=chainlength,initial=0, char="*",style=3,width=wuse)
    
    term=matrix(NA,chainlength,n.models)
    ginverse=ginv[[j]]
    
    for(i in 1:chainlength){   
      cc=post.draw[[j]]()
      psi=ginverse(cc)
      if(i == 101){
        for(k in 1:n.models){
          if(abs(max(detJtest[,k]) - min(detJtest[,k])) < .Machine$double.eps ^ 0.5){samedet[k]=TRUE}
        }
      }
      for(k in 1:n.models){
        gk=g[[k]]
        like=likelihood[[k]]
        prior=param.prior[[k]]
        p=gk(psi)
        if(i<=100) detJtest[i,k] = log(abs(det(attr(adiff(gk, psi), "gradient"))))
        if(samedet[k]==TRUE) detJ=detJtest[1,k] else detJ = log(abs(det(attr(adiff(gk, psi), "gradient"))))
        term[i,k]=like(p)+prior(p)+detJ+model.prior[k]
      }
      term[i,]=term[i,]-max(term[i,])
      term[i,]=exp(term[i,])/sum(exp(term[i,]))
      if(any(is.na(term[i,]))){junk=psi;break}
      utils::setTxtProgressBar(pb,value=i)
    }
    close(pb)
    TM[j,]=apply(term,2,mean)
  }
  ev=eigen(t(TM))
  pi.us=ev$vector[,which(abs(ev$values-1)<1e-8)]
  pi=pi.us/sum(pi.us)
  
  BF=pi/pi[1]*model.prior[1]/model.prior
  
  return(list(TM = TM, prb=pi,BF = BF))
}
