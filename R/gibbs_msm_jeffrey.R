#' Runs a gibbs sampler chain with Jeffreys prior on N to estimate our model parameters
#' 
#' @description Runs a gibbs chain with a jeffrey prior with listing data input 
#' 
#' @param seed a seed
#' @param N number of iterations for the gibbs sampler, default is 10,000
#' @param inclusionProbPriors priors for the inclusion probabilities of participating in each listing. should be present in a configureation file. 
#' @param gibbs.list list containing all the location listing data
#' @param parnames vector containing parameter names to be estimated 
#' 
#' @return a numeric matrix
#' 
#' @importFrom BiasedUrn rFNCHypergeo
#' @importFrom stats rnbinom rbeta rhyper
#' 
#' @export



#note that this gibbs sampler is very specific to our analysis. A fully general 
#gibbs function is beyond me at this point, and if such a thing is possible
#I'm sure that it has been done.  Thus the parameters are actually hardcoded,
#but my hope is that the documentation is good enough that someone could 
#easily adapt the model to their own purposes. 

gibbs_msm_jeffrey = function(seed, N = 10000, inclusionProbPriors, gibbs.list, parnames) {
  
  P.PP = gibbs.list[[1]][1]
  N.srv.PP = gibbs.list[[1]][2]
  N.uid.PP = gibbs.list[[1]][3]
  N.srv.uid.PP = gibbs.list[[1]][4]
  r.PP = N.srv.PP + N.uid.PP - N.srv.uid.PP
  
  P.Nh = gibbs.list[[2]][1]
  N.srv.Nh = gibbs.list[[2]][2]
  N.uid.Nh = gibbs.list[[2]][3]
  N.rnb.Nh = gibbs.list[[2]][4]
  N.srv.uid.Nh = gibbs.list[[2]][5]
  N.srv.rnb.Nh = gibbs.list[[2]][6]
  
  P.ME = gibbs.list[[3]][1]
  N.srv.ME = gibbs.list[[3]][2]
  N.uid.ME = gibbs.list[[3]][3]
  N.srv.uid.ME = gibbs.list[[3]][4]
  N.srv.cpn.ME = gibbs.list[[3]][5]
  
  P.MM = gibbs.list[[4]][1]
  N.srv.MM = gibbs.list[[4]][2]
  N.uid.MM = gibbs.list[[4]][3]
  N.srv.uid.MM = gibbs.list[[4]][4]
  N.srv.cpn.MM = gibbs.list[[4]][5]
  
  N.cpn.Co = gibbs.list[[5]]
  
  probpars=parnames[which(startsWith(parnames,"p"))]
  
  print(seed)
  M=matrix(0,N,length(parnames)) ## stores MCMC samples 
  colnames(M)=parnames
  
  ### constants for beta posteriors ###
  #inclusionProbPriors is meant to be specified in a configuration
  #file created by the user at the start of the analysis.
  a.srv=a.uid=a.rnb=a.cpn=inclusionProbPriors
  b.srv=b.uid=b.rnb=b.cpn=inclusionProbPriors
  
  #set the inclusion probability parameters arbitrarily equal
  #to 0.1
  for(par in probpars) assign(par,0.1)
  
  #initial arbitrary esimate of the population of interest at each location
  #by setting each N.city equal to the number who participated
  #in the survey divided by the probability of participating in that survey
  
  N.Nh=round(N.srv.Nh/p.srv.Nh)
  N.ME=round(N.srv.ME/p.srv.ME)
  N.MM=round(N.srv.MM/p.srv.MM)
  N.cpn.ME=max(N.srv.cpn.ME,round(N.cpn.Co*N.ME/(N.ME+N.MM)))
  N.cpn.MM=N.cpn.Co-N.cpn.ME

  
  set.seed(seed)
  
  for (i in 1:N) {
    
    #### updating Piggs Peak ####
    x=(1-p.srv.PP)*(1-p.uid.PP)
    N.PP=r.PP+rnbinom(1,r.PP,1-x)
    p.srv.PP=rbeta(1,a.srv+N.srv.PP,b.srv+N.PP-N.srv.PP)
    p.uid.PP=rbeta(1,a.uid+N.uid.PP,b.uid+N.PP-N.uid.PP)
    
    ### updating Nhlangano ####
    r.Nh=N.srv.Nh+N.uid.Nh+N.rnb.Nh-N.srv.uid.Nh-N.srv.rnb.Nh - 
      rhyper(1,N.uid.Nh-N.srv.uid.Nh,
             N.Nh-N.srv.Nh-N.uid.Nh+N.srv.uid.Nh,
             N.rnb.Nh-N.srv.rnb.Nh)
    x=(1-p.srv.Nh)*(1-p.uid.Nh)*(1-p.rnb.Nh)
    N.Nh=r.Nh+rnbinom(1,r.Nh,1-x)
    p.srv.Nh=rbeta(1,a.srv+N.srv.Nh,b.srv+N.Nh-N.srv.Nh)
    p.uid.Nh=rbeta(1,a.uid+N.uid.Nh,b.uid+N.Nh-N.uid.Nh)
    p.rnb.Nh=rbeta(1,a.rnb+N.rnb.Nh,b.rnb+N.Nh-N.rnb.Nh)
    
    ### updating Mbabane/Ezulwini ####
    ov.ME=rhyper(1,N.uid.ME-N.srv.uid.ME,N.ME-N.srv.ME-N.uid.ME+N.srv.uid.ME,
                 N.cpn.ME-N.srv.cpn.ME)
    r.ME=N.srv.ME+N.uid.ME+N.cpn.ME-N.srv.uid.ME-N.srv.cpn.ME - ov.ME
    x=(1-p.srv.ME)*(1-p.uid.ME)*(1-p.cpn.ME)
    N.ME=r.ME+rnbinom(1,r.ME,1-x)
    p.srv.ME=rbeta(1,a.srv+N.srv.ME,b.srv+N.ME-N.srv.ME)
    p.uid.ME=rbeta(1,a.uid+N.uid.ME,b.uid+N.ME-N.uid.ME)
    p.cpn.ME=rbeta(1,a.cpn+N.cpn.ME,b.cpn+N.ME-N.cpn.ME)
    
    #### updating Manzini/Matsapha ####
    ov.MM=rhyper(1,N.uid.MM-N.srv.uid.MM,N.MM-N.srv.MM-N.uid.MM+N.srv.uid.MM,N.cpn.MM-N.srv.cpn.MM)
    r.MM=N.srv.MM+N.uid.MM+N.cpn.MM-N.srv.uid.MM-N.srv.cpn.MM - ov.MM
    x=(1-p.srv.MM)*(1-p.uid.MM)*(1-p.cpn.MM)
    N.MM=r.MM+rnbinom(1,r.MM,1-x)
    p.srv.MM=rbeta(1,a.srv+N.srv.MM,b.srv+N.MM-N.srv.MM)
    p.uid.MM=rbeta(1,a.uid+N.uid.MM,b.uid+N.MM-N.uid.MM)
    p.cpn.MM=rbeta(1,a.cpn+N.cpn.MM,b.cpn+N.MM-N.cpn.MM)
    
    #### updating coupon distribution in the corridor ####
    N.cpn.ME=N.srv.cpn.ME+ov.ME+rFNCHypergeo(1,N.ME-N.srv.ME-N.uid.ME+N.srv.uid.ME,
                                             N.MM-N.srv.MM-N.uid.MM+N.srv.uid.MM,
                                             N.cpn.Co-N.srv.cpn.ME-ov.ME-N.srv.cpn.MM-ov.MM,
                                             (p.cpn.ME/(1-p.cpn.ME))/(p.cpn.MM/(1-p.cpn.MM)))
    N.cpn.MM=N.cpn.Co-N.cpn.ME
    
    #even though these aren't used in the simulation we want to 
    #see how well this model can recover them
    
    phi.PP = N.PP/P.PP
    
    phi.Nh = N.Nh/P.Nh
    
    phi.ME = N.ME/P.ME
    
    phi.MM = N.MM/P.MM
    
    for(par in parnames) M[i,par]=get(par)
    
    if((i %% 500)==0) print(i)    
    #print(i)
    
  }
  

  #we burn the first half of the gibbs iterations
  M.burn = M[(nrow(M)/2):nrow(M),]
  
  return(M.burn)
}