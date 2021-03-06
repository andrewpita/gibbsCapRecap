#' Runs a gibbs sampler chain with Jeffreys prior on N to estimate our model parameters
#' for the FSW data
#' 
#' @description Runs a gibbs chain with a jeffrey prior with listing data input
#' 
#' @param seed a seed
#' @param N number of iterations for the gibbs sampler. Default is 10,000
#' @param inclusionProbPriors priors for the inclusion probabilities of participating in each listing. should be present in a configureation file. 
#' @param constraint Value that represents the lower bound "known" number of flas participants in Mbabane
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

gibbs_fsw_jeffrey = function(seed, N = 10000, inclusionProbPriors,
                             constraint,gibbs.list, parnames) {
  
  P.PP = gibbs.list[[1]][1]
  N.srv.PP = gibbs.list[[1]][2]
  N.uid.PP = gibbs.list[[1]][3]
  N.srv.uid.PP = gibbs.list[[1]][4]
  r.PP = N.srv.PP + N.uid.PP - N.srv.uid.PP
  
  P.Lv = gibbs.list[[2]][1]
  N.srv.Lv = gibbs.list[[2]][2]
  N.uid.Lv = gibbs.list[[2]][3]
  N.srv.uid.Lv = gibbs.list[[2]][4]
  r.Lv = N.srv.Lv + N.uid.Lv - N.srv.uid.Lv
  
  P.Nh = gibbs.list[[3]][1]
  N.srv.Nh = gibbs.list[[3]][2]
  N.uid.Nh = gibbs.list[[3]][3]
  N.srv.uid.Nh = gibbs.list[[3]][4]
  r.Nh = N.srv.Nh + N.uid.Nh - N.srv.uid.Nh
  
  P.ME = gibbs.list[[4]][1]
  N.srv.ME = gibbs.list[[4]][2]
  N.uid.ME = gibbs.list[[4]][3]
  N.srv.uid.ME = gibbs.list[[4]][4]
  N.srv.flas.ME = gibbs.list[[4]][5]
  
  P.MM = gibbs.list[[5]][1]
  N.srv.MM = gibbs.list[[5]][2]
  N.uid.MM = gibbs.list[[5]][3]
  N.srv.uid.MM = gibbs.list[[5]][4]
  N.srv.flas.MM = gibbs.list[[5]][5]
  
  N.flas.Co = gibbs.list[[6]]
  
  probpars=parnames[which(startsWith(parnames,"p"))]
  
  print(seed)
  M=matrix(0,N,length(parnames)) ## stores MCMC samples 
  colnames(M)=parnames
  
  ### constants for beta posteriors ###
  #inclusionProbPriors is meant to be specified in a configuration
  #file created by the user at the start of the analysis.
  a.srv=a.uid=a.flas=inclusionProbPriors
  b.srv=b.uid=b.flas=inclusionProbPriors
  
  #set the inclusion probability parameters arbitrarily equal
  #to 0.1
  for(par in probpars) assign(par,0.1)
  
  #initial arbitrary esimate of the population of interest at each location
  #by setting each N.city equal to the number who participated
  #in the survey divided by the probability of participating in that survey
  
  N.Lv = round(N.srv.Lv/p.srv.Lv)
  N.PP = round(N.srv.PP/p.srv.PP)
  N.Nh=round(N.srv.Nh/p.srv.Nh)
  N.ME=round(N.srv.ME/p.srv.ME)
  N.MM=round(N.srv.MM/p.srv.MM)
  N.flas.ME=max(N.srv.flas.ME,round(N.flas.Co*N.ME/(N.ME+N.MM)))
  N.flas.MM=N.flas.Co-N.flas.ME
  
  set.seed(seed)
  
  for (i in 1:N) {
    
    #updating Lavumisa
    x = (1 - p.srv.Lv)*(1 - p.uid.Lv)
    N.Lv = r.Lv + rnbinom(1, r.Lv, 1 - x)
    p.srv.Lv  = rbeta(1, a.srv + N.srv.Lv, b.srv + N.Lv - N.srv.Lv)
    p.uid.Lv = rbeta(1, a.uid + N.uid.Lv, b.srv + N.Lv - N.uid.Lv)
    
    #updating Piggs Peak
    x=(1-p.srv.PP)*(1-p.uid.PP)
    N.PP=r.PP+rnbinom(1,r.PP,1-x)
    p.srv.PP=rbeta(1,a.srv+N.srv.PP,b.srv+N.PP-N.srv.PP)
    p.uid.PP=rbeta(1,a.uid+N.uid.PP,b.uid+N.PP-N.uid.PP)
    
    #updating Nhlangano
    x = (1 - p.srv.Nh)*(1 - p.uid.Nh)
    N.Nh = r.Nh + rnbinom(1, r.Nh, 1 - x)
    p.srv.Nh = rbeta(1, a.srv + N.srv.Nh, b.srv + N.Nh - N.srv.Nh)
    p.uid.Nh = rbeta(1, a.uid + N.uid.Nh, b.uid + N.Nh - N.uid.Nh)
    
    #updating Mbabane/Ezulwini
    ov.ME = rhyper(1, N.uid.ME - N.srv.uid.ME, N.ME - N.srv.ME - N.uid.ME + N.srv.uid.ME,
                   N.flas.ME - N.srv.flas.ME)
    r.ME = N.srv.ME + N.uid.ME + N.flas.ME - N.srv.uid.ME - N.srv.flas.ME - ov.ME
    x = (1 - p.srv.ME)*(1 - p.uid.ME)*(1 - p.flas.ME)
    N.ME = r.ME + rnbinom(1, r.ME, 1-x)
    p.srv.ME = rbeta(1, a.srv + N.srv.ME, b.srv + N.ME - N.srv.ME)
    p.uid.ME = rbeta(1, a.uid + N.uid.ME, b.srv + N.ME - N.uid.ME)
    p.flas.ME = rbeta(1, a.flas + N.flas.ME, b.srv + N.ME - N.flas.ME)
    
    #updating Manzini/Matsapha
    ov.MM = rhyper(1, N.uid.MM - N.srv.uid.MM, N.MM - N.srv.MM - N.uid.MM + N.srv.uid.MM,
                   N.flas.MM - N.srv.flas.MM)
    r.MM = N.srv.MM + N.uid.MM + N.flas.MM - N.srv.uid.MM - N.srv.flas.MM - ov.MM
    x = (1 - p.srv.MM)*(1 - p.uid.MM)*(1 - p.flas.MM)
    N.MM = r.MM + rnbinom(1, r.MM, 1-x)
    p.srv.MM = rbeta(1, a.srv + N.srv.MM, b.srv + N.MM - N.srv.MM)
    p.uid.MM = rbeta(1, a.uid + N.uid.MM, b.srv + N.MM - N.uid.MM)
    p.flas.MM = rbeta(1, a.flas + N.flas.MM, b.srv + N.MM - N.flas.MM)
    
    #updating corridor
    
    temp = 0
    
    while (temp < constraint) {
      
      N.flas.ME = N.srv.flas.ME + ov.ME + rFNCHypergeo(1, N.ME - N.srv.ME - N.uid.ME + N.srv.uid.ME,
                                                       N.MM - N.srv.MM - N.uid.MM + N.srv.uid.MM,
                                                       N.flas.Co - N.srv.flas.ME - ov.ME - N.srv.flas.MM - ov.MM,
                                                       ((p.flas.ME)/(1-p.flas.ME))/((p.flas.MM)/(1-p.flas.MM)))
      temp = N.flas.ME
      
    }
    
    N.flas.MM = N.flas.Co - N.flas.ME
    
    phi.Lv = N.Lv/P.Lv
    
    phi.PP = N.PP/P.PP
    
    phi.Nh = N.Nh/P.Nh
    
    phi.ME = N.ME/P.ME
    
    phi.MM = N.MM/P.MM
    
    #store this iterations values for all parameters
    for(par in parnames) M[i,par]=get(par)
    
    if((i %% 500)==0) print(i)
    
  }
  
  #we burn the first half of the gibbs iterations
  M.burn = M[(nrow(M)/2):nrow(M),]
  
  return(M.burn)
  
}