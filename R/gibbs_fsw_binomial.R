#' Runs a gibbs sampler chain with a Binomial prior on N to estimate our model parameters
#' for FSW data
#' 
#' @description Requires that parameters be declared globally, and note they are 
#' @description hardcoded for this analysis
#' 
#' @param seed a seed
#' @param N number of iterations for the gibbs sampler. Default 10,000
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




gibbs_fsw_binomial = function(seed, N = 10000, atune = 0.25,
                              btune = 1, inclusionProbPriors) {
  
  print(seed)
  
  M=matrix(0,N,length(parnames)) ## stores MCMC samples 
  colnames(M)=parnames
  
  
  a.srv = a.uid = a.flas = inclusionProbPriors
  b.srv = b.uid = b.flas = inclusionProbPriors
  
  ### initializing the chain ###
  #initial fsw population estimates are the number of people
  #at each location who took the survey, divided by the initial 
  #probability of taking the survey
  for(par in probpars) assign(par,0.1)
  
  N.Lv = round(P.Lv*phi.Lv)
  N.PP = round(P.PP*phi.PP)
  N.Nh=round(P.Nh*phi.Nh)
  N.ME=round(P.ME*phi.ME)
  N.MM=round(P.MM*phi.MM)
  N.flas.ME=max(N.srv.flas.ME,round(N.flas.Co*N.ME/(N.ME+N.MM)))
  N.flas.MM=N.flas.Co-N.flas.ME
  a.phi=b.phi=1.5
  
  set.seed(seed)
  
  for (i in 1:N) {
    
    #### updating Lavumisa ####
    x=phi.Lv*(1-p.srv.Lv)*(1-p.uid.Lv)
    N.Lv=r.Lv+rbinom(1,P.Lv-r.Lv,x/(1+x))
    p.srv.Lv=rbeta(1,a.srv+N.srv.Lv,b.srv+N.Lv-N.srv.Lv)
    p.uid.Lv=rbeta(1,a.uid+N.uid.Lv,b.uid+N.Lv-N.uid.Lv)
    phi.Lv=rbeta(1,a.phi+N.Lv,b.phi+P.Lv-N.Lv)
    
    #### updating Piggs Peak ####
    x=phi.PP*(1-p.srv.PP)*(1-p.uid.PP)
    N.PP=r.PP+rbinom(1,P.PP-r.PP,x/(1+x))
    p.srv.PP=rbeta(1,a.srv+N.srv.PP,b.srv+N.PP-N.srv.PP)
    p.uid.PP=rbeta(1,a.uid+N.uid.PP,b.uid+N.PP-N.uid.PP)
    phi.PP=rbeta(1,a.phi+N.PP,b.phi+P.PP-N.PP)
    
    #### updating Nhlangano ####
    x=phi.Nh*(1-p.srv.Nh)*(1-p.uid.Nh)
    N.Nh=r.Nh+rbinom(1,P.Nh-r.Nh,x/(1+x))
    p.srv.Nh=rbeta(1,a.srv+N.srv.Nh,b.srv+N.Nh-N.srv.Nh)
    p.uid.Nh=rbeta(1,a.uid+N.uid.Nh,b.uid+N.Nh-N.uid.Nh)
    phi.Nh=rbeta(1,a.phi+N.Nh,b.phi+P.Nh-N.Nh)
    
    ### updating Mbabane/Ezulwini ####
    ov.ME=rhyper(1,N.uid.ME-N.srv.uid.ME,N.ME-N.srv.ME-N.uid.ME+N.srv.uid.ME,
                 N.flas.ME-N.srv.flas.ME)
    r.ME=N.srv.ME+N.uid.ME+N.flas.ME-N.srv.uid.ME-N.srv.flas.ME - ov.ME
    x=phi.ME*(1-p.srv.ME)*(1-p.uid.ME)*(1-p.flas.ME)
    N.ME=r.ME+rbinom(1,P.ME-r.ME,x/(1+x))
    p.srv.ME=rbeta(1,a.srv+N.srv.ME,b.srv+N.ME-N.srv.ME)
    p.uid.ME=rbeta(1,a.uid+N.uid.ME,b.uid+N.ME-N.uid.ME)
    p.flas.ME=rbeta(1,a.flas+N.flas.ME,b.flas+N.ME-N.flas.ME)
    phi.ME=rbeta(1,a.phi+N.ME,b.phi+P.ME-N.ME)
    
    
    #### updating Manzini/Matsapha ####
    ov.MM=rhyper(1,N.uid.MM-N.srv.uid.MM,N.MM-N.srv.MM-N.uid.MM+N.srv.uid.MM,
                 N.flas.MM-N.srv.flas.MM)
    r.MM=N.srv.MM+N.uid.MM+N.flas.MM-N.srv.uid.MM-N.srv.flas.MM - ov.MM
    x=phi.MM*(1-p.srv.MM)*(1-p.uid.MM)*(1-p.flas.MM)
    N.MM=r.MM+rbinom(1,P.MM-r.MM,x/(1+x))
    p.srv.MM=rbeta(1,a.srv+N.srv.MM,b.srv+N.MM-N.srv.MM)
    p.uid.MM=rbeta(1,a.uid+N.uid.MM,b.uid+N.MM-N.uid.MM)
    p.flas.MM=rbeta(1,a.flas+N.flas.MM,b.flas+N.MM-N.flas.MM)
    phi.MM=rbeta(1,a.phi+N.MM,b.phi+P.MM-N.MM)
    
    #updating Corridor
    
    temp = 0
    
    while (temp < 70) {
      
      N.flas.ME = N.srv.flas.ME + ov.ME + 
        
        rFNCHypergeo(1, N.ME - N.srv.ME - N.uid.ME + N.srv.uid.ME,
                     N.MM - N.srv.MM - N.uid.MM + N.srv.uid.MM,
                     N.flas.Co - N.srv.flas.ME - ov.ME - N.srv.flas.MM - ov.MM,
                     (p.flas.ME)/(1-p.flas.ME)/(p.flas.MM)/(1-p.flas.MM))
      temp = N.flas.ME
      
    }
    
    N.flas.MM = N.flas.Co - N.flas.ME
    
    #### Metropolis update for a.phi and b.phi ####
    abphi=betahyperMH(a.phi,b.phi,atune,btune,c(phi.Lv,phi.PP,phi.Nh,phi.ME,phi.MM))
    #print(abphi)
    a.phi=abphi[1]
    b.phi=abphi[2]
    
    for(par in parnames) M[i,par]=get(par)
    
    if((i %% 1)==200) print(i)
    
    
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