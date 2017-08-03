
#' Methods for the rj Class
#' 
#' An object of class \code{rj} is returned from the functions 
#' \code{\link{rjmcmcpost}} or \code{\link{defaultpost}}. The following methods
#' can be applied to an object of this class. See Details for more information.
#' 
#' The \code{print} method prints the point estimates obtained from the
#' algorithm, including the transition matrix, posterior model probabilities and
#' Bayes factors.
#' 
#' The \code{probplot} function plots how the estimates of the posterior probabilities
#' changed as the algorithm progressed, illustrating convergence.
#' 
#' The \code{densities} function returns the likelihood, prior and posterior
#' density for each model at each iteration of the algorithm in a list of length K, where each element is a 3K x
#' \code{chainlength} matrix. 
#' 
#' The \code{psidraws} function returns the universal
#' parameter vector psi at each iteration of the algorithm in a matrix.  
#' 
#' The \code{detJ} function returns the logarithm of the determinant of the Jacobian
#' matrix for each model, in a vector of length K.
#' 
#' @param rjobj,x An object of class \code{rj}.
#' @param ... Any further arguments to \code{print}.
#' @param legend,col,ylim,lwd,lty Optional graphical parameters to the \code{probplot} function.
#' 
#' @importFrom graphics lines plot
#' @name rjmethods
NULL

#' @rdname rjmethods
#' @export
print.rj = function(x, ...){
  if(any(names(x)=="result")){ print(x$result) } else { warning("Object does not contain an element named 'result'.") }
}

#' @rdname rjmethods
#' @export
probplot = function(rjobj, legend=TRUE, col="maroon4", ylim=c(0,1), lwd=2, lty=c(1,1,1)){
  numplots = dim(rjobj$progress$prb)[2]
  if(length(col)==1){ col = c(col, 2:numplots+2) }
  plot(seq(rjobj$meta$TM.thin,rjobj$meta$chainl,by=rjobj$meta$TM.thin), rjobj$progress$prb[,1], type="l", lty=lty, lwd=lwd, col=col[1], ylim=ylim, main="Estimated Posterior Model Probabilities", rjobjlab="Number of Iterations", ylab="Estimated Posterior Probability")
  for(i in 2:numplots){
    lines(seq(rjobj$meta$TM.thin,rjobj$meta$chainl,by=rjobj$meta$TM.thin), rjobj$progress$prb[,i], lwd=lwd, col=col[i])
  }
  if(legend){ legend("topright", legend=paste("p(Model", 1:numplots, "| y)"), lwd=2, col=col) }
}

#' @rdname rjmethods
#' @export
densities = function(rjobj){
  if(any(names(rjobj)=="densities")){ return(rjobj$densities) } else { warning("Object does not contain an element named 'densities'.") }
}

#' @rdname rjmethods
#' @export
psidraws = function(rjobj){
  if(any(names(rjobj)=="psidraws")){ return(rjobj$psidraws) } else { warning("Object does not contain an element named 'psidraws'.") }
}

#' @rdname rjmethods
#' @export
detJ = function(rjobj){
  if(any(names(rjobj)=="detJ")){ return(rjobj$detJ) } else { warning("Object does not contain an element named 'detJ'.") }
}

rj = function(rjobj){
  class(rjobj) = c(class(rjobj), "rj")
  return(rjobj)
}
