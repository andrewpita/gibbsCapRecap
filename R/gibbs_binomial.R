#' Runs a gibbs sampler chain with a Binomial prior on N to estimate our model parameters
#' 
#' @description Requires that parameters be declared globally, and note they are 
#' @description hardcoded for this analysis
#' 
#' @param seed a seed
#' @param N number of iterations for the gibbs sampler
#' @param inclusionProbPriors priors for the inclusion probabilities of participating in each listing. should be present in a configureation file. 
#' @param atune Tuning parameter for the metropolis-hastings update of hyperparameters on phis
#' @param btune Other tuning parameter for the metropolis-hastings update
#' @param x log of a.phi/(a.phi+ b.phi)
#' @param y log of (a.phi +b.phi)
#' @param phivec vector of current values of phi
#' 
#' @return a numeric matrix
#' 
#' @importFrom BiasedUrn rFNCHypergeo
#' @importFrom stats rnbinom rbeta rhyper runif rbinom dbeta
#' @importFrom truncnorm rtruncnorm
#' 
#' @export



gibbs_binomial = function(seed, N = 10000, atune = 0.25,
                          btune = 1, inclusionProbPriors) {
  
  print(seed)
  M=matrix(0,N,length(parnames)) ## stores MCMC samples 
  colnames(M)=parnames
  
  ### constants for beta posteriors ###
  a.srv=a.uid=a.rnb=a.cpn=inclusionProbPriors
  b.srv=b.uid=b.rnb=b.cpn=inclusionProbPriors
  
  ### initializing the chain ###
  #arbitrary initial estimates fo inclusion probabilities and 
  #target population proportions is 0.1
  for(par in probpars) assign(par,0.1)
  
  #initial target population estimates are the number of people
  #at each location who could have taken the survey, 
  #multiplied by the initial estimate of the proportion 
  #of the population that is the target population
  N.Nh=round(P.Nh*phi.Nh)
  N.ME=round(P.ME*phi.ME)
  N.MM=round(P.MM*phi.MM)
  N.cpn.ME=max(N.srv.cpn.ME,round(N.cpn.Co*N.ME/(N.ME+N.MM)))
  N.cpn.MM=N.cpn.Co-N.cpn.ME
  a.phi=b.phi=1.5
  
  set.seed(seed)
  
  for (i in 1:N) {
    
    #### updating Piggs Peak ####
    x=phi.PP*(1-p.srv.PP)*(1-p.uid.PP)
    N.PP=r.PP+rbinom(1,P.PP-r.PP,x/(1+x))
    p.srv.PP=rbeta(1,a.srv+N.srv.PP,b.srv+N.PP-N.srv.PP)
    p.uid.PP=rbeta(1,a.uid+N.uid.PP,b.uid+N.PP-N.uid.PP)
    phi.PP=rbeta(1,a.phi+N.PP,b.phi+P.PP-N.PP)
    
    ### updating Nhlangano ####
    r.Nh=N.srv.Nh+N.uid.Nh+N.rnb.Nh-N.srv.uid.Nh-N.srv.rnb.Nh - 
      rhyper(1,N.uid.Nh-N.srv.uid.Nh,N.Nh-N.srv.Nh-N.uid.Nh+N.srv.uid.Nh,N.rnb.Nh-N.srv.rnb.Nh)
    x=phi.Nh*(1-p.srv.Nh)*(1-p.uid.Nh)*(1-p.rnb.Nh)
    N.Nh=r.Nh+rbinom(1,P.Nh-r.Nh,x/(1+x))
    p.srv.Nh=rbeta(1,a.srv+N.srv.Nh,b.srv+N.Nh-N.srv.Nh)
    p.uid.Nh=rbeta(1,a.uid+N.uid.Nh,b.uid+N.Nh-N.uid.Nh)
    p.rnb.Nh=rbeta(1,a.rnb+N.rnb.Nh,b.rnb+N.Nh-N.rnb.Nh)
    phi.Nh=rbeta(1,a.phi+N.Nh,b.phi+P.Nh-N.Nh)
    
    ### updating Mbabane/Ezulwini ####
    ov.ME=rhyper(1,N.uid.ME-N.srv.uid.ME,N.ME-N.srv.ME-N.uid.ME+N.srv.uid.ME,N.cpn.ME-N.srv.cpn.ME)
    r.ME=N.srv.ME+N.uid.ME+N.cpn.ME-N.srv.uid.ME-N.srv.cpn.ME - ov.ME
    x=phi.ME*(1-p.srv.ME)*(1-p.uid.ME)*(1-p.cpn.ME)
    N.ME=r.ME+rbinom(1,P.ME-r.ME,x/(1+x))
    p.srv.ME=rbeta(1,a.srv+N.srv.ME,b.srv+N.ME-N.srv.ME)
    p.uid.ME=rbeta(1,a.uid+N.uid.ME,b.uid+N.ME-N.uid.ME)
    p.cpn.ME=rbeta(1,a.cpn+N.cpn.ME,b.cpn+N.ME-N.cpn.ME)
    phi.ME=rbeta(1,a.phi+N.ME,b.phi+P.ME-N.ME)
    
    #### updating Manzini/Matsapha ####
    ov.MM=rhyper(1,N.uid.MM-N.srv.uid.MM,N.MM-N.srv.MM-N.uid.MM+N.srv.uid.MM,N.cpn.MM-N.srv.cpn.MM)
    r.MM=N.srv.MM+N.uid.MM+N.cpn.MM-N.srv.uid.MM-N.srv.cpn.MM - ov.MM
    x=phi.MM*(1-p.srv.MM)*(1-p.uid.MM)*(1-p.cpn.MM)
    N.MM=r.MM+rbinom(1,P.MM-r.MM,x/(1+x))
    p.srv.MM=rbeta(1,a.srv+N.srv.MM,b.srv+N.MM-N.srv.MM)
    p.uid.MM=rbeta(1,a.uid+N.uid.MM,b.uid+N.MM-N.uid.MM)
    p.cpn.MM=rbeta(1,a.cpn+N.cpn.MM,b.cpn+N.MM-N.cpn.MM)
    phi.MM=rbeta(1,a.phi+N.MM,b.phi+P.MM-N.MM)
    
    #### updating coupon distribution in the corridor ####
    N.cpn.ME=N.srv.cpn.ME+ov.ME+rFNCHypergeo(1,N.ME-N.srv.ME-N.uid.ME+N.srv.uid.ME,
                                             N.MM-N.srv.MM-N.uid.MM+N.srv.uid.MM,N.cpn.Co-N.srv.cpn.ME-ov.ME-N.srv.cpn.MM-ov.MM,
                                             (p.cpn.ME/(1-p.cpn.ME))/(p.cpn.MM/(1-p.cpn.MM)))
    N.cpn.MM=N.cpn.Co-N.cpn.ME
    
    #### Metropolis update for a.phi and b.phi ####
    abphi=betahyperMH(a.phi,b.phi,atune,btune,c(phi.PP,phi.Nh,phi.ME,phi.MM))
    #print(abphi)
    a.phi=abphi[1]
    b.phi=abphi[2]
    
    for(par in parnames) M[i,par]=get(par)
    
    if((i %% 1)==200) print(i)    
    #print(i)
    
  }
  
  #add an overal posterior mean to matrix
  Mphi=as.matrix(M[,"a.phi"]/(M[,"a.phi"]+M[,"b.phi"]))
  colnames(Mphi)="phi"
  M = cbind(M,Mphi)
  
  #burn first half of 
  M.burn = M[(nrow(M)/2):nrow(M),]
  
}


## posterior likelihood ##
postlike=function(x,y,phivec){
  sum(sapply(phivec,dbeta,exp(x+y),(1-exp(x))*exp(y),log=TRUE))+x
}

### MH step ###
betahyperMH=function(a,b,atune,btune,phivec){
  #print(phivec)
  x=log(a/(a+b))
  y=log(a+b)
  
  xnew=rtruncnorm(1,a=-y,b=log(1-exp(-y)),x,atune)
  rx=postlike(xnew,y,phivec)-postlike(x,y,phivec)
  #print(xnew)
  
  ux=runif(1,0,1)
  if(rx>log(ux)) x=xnew
  
  ynew=rtruncnorm(1,a=max(-x,-log(1-exp(x))),b=Inf,y,btune)
  ry=postlike(x,ynew,phivec)-postlike(x,y,phivec)
  #print(ynew)
  
  uy=runif(1,0,1)
  if(ry>log(uy)) y=ynew
  
  #print(c(x,y))
  c(exp(x+y),(1-exp(x))*exp(y))
  
}