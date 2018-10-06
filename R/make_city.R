#' Get the values for parameters for each location for use in simulation.
#' 
#' @description To simulate data we need values for parameters at each location 
#' 
#' @param string The output of the Gibbs Chain contains parameter estimates labeled by location. This string corresponds to that location. 
#' @param df a named vector or df that contains the parameter values
#' 
#' @return a numeric vector
#' 


#in our analysis we  specifically used the median of the posterior distribution approximated using the 
#gibbs sampler.  Thus for df we use a named vector of the median estimates for each of the parameters.
#the last two characters of the parameter names correspond to a specific location

make_city = function(string, df) {
  
  indeces = grepl(paste(string,"$", sep = ""), names(df))
  
  city = df[indeces]
  names(city) = NULL
  return(city)
  
}



