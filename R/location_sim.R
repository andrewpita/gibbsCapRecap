#' Simulates draws from multinomial model of listings for 2 or 3 listings
#' 
#' @description works for 2 or 3 listings, probs are the inclusion probabilites  
#' @description for each listing
#' 
#' @param nListings number of listings, 2 or 3
#' @param size total size drawn
#' @param probs vector of inclusion probabilities for each listing
#' 
#' @return a list of draws
#' 
#' @importFrom stats rmultinom
#' 
#' @export


location_sim = function(nListings, size, probs) {
  
  if (nListings == 2) {
    
    prob1 = probs[1]
    prob2 = probs[2]
    
    samp = rmultinom(1, size = size, prob = 
                       c(prob1 * (1 - prob2),#[1] in first listing
                         prob1 * prob2,#[2] in both
                         (1 - prob1) * prob2,#[3] in second only
                         (1- prob1) * (1 - prob2)#[4] in neither
                       ))
    N.both = samp[2]
    N.1 = samp[1] + N.both
    N.2 = samp[2] + samp[3]
    r = N.1 + N.2 - N.both #in at least one
    r.list = list(N.1, N.2, N.both, r)
    return(r.list)
    
  } else if (nListings == 3) {
    
    prob1 = probs[1]
    prob2 = probs[2]
    prob3 = probs[3]
    
    samp = rmultinom(1, size = size, prob = 
                       c(prob1 * (1 - prob2) * (1 - prob3),#[1] in only first
                         prob1 * prob2 * (1- prob3),#[2] in 1 and 2 only
                         prob1 * prob2 * prob3,#[3] in all three
                         (1- prob1) * prob2 * (1- prob3),#[4] in only 2
                         (1 - prob1) * prob2 * prob3,#[5] in 2 and 3 only
                         (1 - prob1) * (1 - prob2) * prob3,#[6] in only 3
                         prob1 * (1 - prob2) * prob3,#[7] in 1 and 3 only
                         (1 - prob1) * (1 -  prob2) * (1 -  prob3)#[8] in none
                       ))
    
    N.1.2 = samp[2] + samp[3]
    N.1.3 = samp[3] + samp[7]
    N.1 = samp[1] + samp[2] + samp[3] + samp[7]
    N.2 = samp[2] + samp[3] + samp[4] + samp[5]
    N.3 = samp[3] + samp[5] + samp[6] + samp[7]
    r = as.numeric(size - samp[8])
    r.list = list(N.1, N.2, N.3, N.1.2, N.1.3, r)
    return(r.list)
    
  } else {
    
    print("nListings must be equal to 2 or 3")
  }
  
}