#' Simulates data for msm analysis and uses gibbs_msm_binomial
#' to evaluate how often the chain draws cover the true values
#' 
#' @description Requires that parameters be declared globally, and note they are 
#' @description hardcoded for this analysis
#' 
#' @param seed a seed
#' @param N number of iterations for the gibbs sampler, default is 10,000
#' @param inclusionProbPriors priors for the inclusion probabilities of participating in each listing. should be present in a configureation file. 
#' @param atune tuning parameter for metropolis-hastings step
#' @param btune tuning parameter for metropolis-hastings step
#' 
#' @return a list that contains one iterations gibbs chain, a summary, and whether or not each parameter was recovered
#' 
#' @importFrom stats rbinom rmultinom
#' 
#' @export

sim_msm_binomial = function(seed, N = 10000, inclusionProbPriors,
                            atune = 0.25, btune = 1) {
  
  #simulate data for Pigg's Peak
  
  true.N.PP = rbinom(size = P.PP, n = 1, prob = true.phi.PP)
  
  PP.sample = rmultinom(1, size = true.N.PP, prob = 
                          c(true.p.srv.PP*(1 - true.p.uid.PP), #in survey not in uid
                            true.p.srv.PP*true.p.uid.PP, #in both
                            true.p.uid.PP*(1 - true.p.srv.PP), #in uid not in survey
                            (1 - true.p.srv.PP)*(1 - true.p.srv.PP) )) #in neither
  
  #store Pigg's Peak simulated data
  N.srv.uid.PP = PP.sample[2]
  N.srv.PP = PP.sample[1] + N.srv.uid.PP
  N.uid.PP = PP.sample[3] + N.srv.uid.PP
  r.PP = N.srv.PP + N.uid.PP - N.srv.uid.PP
  
  #simulate data for Nhlangano 
  
  true.N.Nh = rbinom(size = P.Nh, n = 1, prob = true.phi.Nh)
  
  Nh.sample = rmultinom(1, size = true.N.Nh, prob = 
                          c(true.p.srv.Nh * (1 - true.p.uid.Nh) * (1 - true.p.rnb.Nh),#[1] in survey only
                            true.p.srv.Nh * true.p.uid.Nh * (1 - true.p.rnb.Nh), #[2]in survey and uid only
                            true.p.srv.Nh * true.p.uid.Nh * true.p.rnb.Nh), #[3]in all three
                        true.p.uid.Nh * (1- true.p.srv.Nh) * (1 - true.p.rnb.Nh),#[4]in only uid
                        true.p.uid.Nh * (1- true.p.srv.Nh) * true.p.rnb.Nh,#[5]in uid and rainbow only
                        true.p.rnb.Nh * (1 - true.p.srv.Nh) * (1 - true.p.uid.Nh),#[6]in only rainbow
                        true.p.rnb.Nh * true.p.srv.Nh * (1- true.p.uid.Nh),#[7]in rainbow and survey only
                        (1- true.p.srv.Nh) * (1- true.p.uid.Nh) * (1 - true.p.rnb.Nh))#[8]in none
  
  #storing simulated data for Nhlangano
  N.srv.uid.Nh = Nh.sample[2] + Nh.sample[3] 
  N.srv.rnb.Nh = Nh.sample[7] + Nh.sample[3]
  N.srv.Nh = Nh.sample[1] + Nh.sample[2] + Nh.sample[3] + Nh.sample[7]
  N.uid.Nh = Nh.sample[2] + Nh.sample[3] + Nh.sample[4] + Nh.sample[5]
  N.rnb.Nh = Nh.sample[3] + Nh.sample[5] + Nh.sample[6] + Nh.sample[7]
  true.r.Nh = as.numeric(true.N.Nh - Nh.sample[8])
  
  #simulate data for Mbabane/Ezulwini
  
  true.N.ME = rbinom(size = P.ME, n = 1, prob = true.phi.ME)
  
  ME.sample = rmultinom(1, size = true.N.ME, prob = 
                          c(true.p.srv.ME * (1- true.p.uid.ME) * (1 - true.p.cpn.ME),#[1]in only survey
                            true.p.srv.ME * true.p.uid.ME * (1 - true.p.cpn.ME),#[2]in survey and uid only
                            true.p.srv.ME * true.p.uid.ME * true.p.cpn.ME,#[3] in all three
                            true.p.uid.ME * (1 - true.p.srv.ME) * (1 - true.p.cpn.ME),#[4] in only uid
                            true.p.uid.ME * (1 - true.p.srv.ME) * (true.p.cpn.ME),#[5] in uid and coupon only
                            true.p.cpn.ME * (1 - true.p.srv.ME) * (1 - true.p.uid.ME),#[6] in only coupon
                            true.p.cpn.ME * true.p.srv.ME * (1 - true.p.uid.ME),#[7] in coupon and survey only
                            (1 - true.p.srv.ME) * (1- true.p.uid.ME) * (1- true.p.cpn.ME)#[8] in none
                          ))
  
  #storing simulated data for Mbabane/Ezulwini
  true.N.cpn.ME = as.numeric(ME.sample[3] + ME.sample[5] + ME.sample[6] +
                               ME.sample[7])
  
  N.srv.uid.ME = ME.sample[2] + ME.sample[3]
  N.srv.cpn.ME = ME.sample[3] + ME.sample[7]
  N.srv.ME = ME.sample[1] + ME.sample[2] + ME.sample[3] + ME.sample[7]
  N.uid.ME = ME.sample[2] + ME.sample[3] + ME.sample[4] + ME.sample[5]
  true.r.ME = as.numeric(true.N.ME - ME.sample[8])
  
  #simulate data for Manzni/Matsapha
  true.N.MM = rbind(size = P.MM, n = 1, prob = true.phi.MM)
  
  MM.sample = rmultinom(1, size = true.N.MM, prob = 
                          c(true.p.srv.MM * (1 - true.p.uid.MM) * (1 - true.p.cpn.MM),#[1] in only survey
                            true.p.srv.MM * true.p.uid.MM * (1 - true.p.cpn.MM),#[2] in survey and uid
                            true.p.srv.MM * true.p.uid.MM * true.p.cpn.MM,#[3] in all three
                            true.p.uid.MM * (1 - true.p.srv.MM) * (1 - true.p.cpn.MM),#[4] in only uid
                            true.p.uid.MM * (1 - true.p.srv.MM) * true.cpn.MM,#[5] in only uid and coupon
                            true.p.cpn.MM * (1 - true.p.srv.MM) * (1 - true.p.uid.MM),#[6] in only coupon
                            true.p.cpn.MM * true.p.srv.MM * (1 - true.p.uid.MM),#[7] in coupon and survey
                            (1 - true.p.cpn.MM) * (1 - true.p.srv.MM) * (1 - true.p.uid.MM)#[8] in none
                          ))
  
  #store simulated data for Manzini/Matsapha
  
  true.N.cpn.MM = MM.sample[3] + MM.sample[5] + MM.sample[6] + MM.sample[7]
  
  N.srv.uid.MM = MM.sample[2] + MM.sample[3]
  N.srv.cpn.MM = MM.sample[3] + MM.sample[7]
  N.srv.MM = MM.sample[1] + MM.sample[2] + MM.sample[3] + MM.sample[7]
  N.uid.MM = MM.sample[2] + MM.sample[3] + MM.sample[4] + MM.sample[5]
  true.r.MM = as.numeric(true.N.MM - MM.sample[8])
  
  # Corridor: sum of MM and ME
  N.cpn.Co = true.N.cpn.MM + true.N.cpn.ME
  
  #we run just one gibbs chain with 10,000 iterations
  M = gibbs_msm_binomial(seed,N, inclusionProbPriors)                     
  
  #apply the gibbs_summary function to each column of M
  summary = apply(M,2,sim_summary)
  
  #since the true values are declared globally using output
  #from analysis but some need to be simulated and are then
  #declared within the function, we grab them from the parent frame
  
  truth=sapply(paste0("true.",parnames),get,envir=sys.frame(sys.parent(0)), simplify=FALSE) 
  names(truth)=colnames(summary)
  
  #a 1 or 0 value for each parameter that tells whether or 
  #not it was contained in the 95% confidence interval of the 
  #gibbs samples
  cov=matrix((summary["2.5%",]< truth)*(summary["97.5%",]>truth),nrow=1)
  row.names(cov)="coverage"
  
  return.list = list(M, rbind(summary, cov), truth)
  
  return(return.list)
  
}