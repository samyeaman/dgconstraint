#' Estimate the effective proportion of adaptive loci
#'
#' @description \code{estimate_pa} Estimates the effective proportion of
#' adaptive loci (P_a) using likelihood-based tests. Setting na.rm = T
#' assumes that missing values are true negatives.
#'
#' @param input numeric array
#' @param ndigits numeric
#' @param show.plot boolean
#'
#' @return length-one numeric
#'
#' @examples
#'
#'
#' @export


estimate_pa <- function(input, ndigits = 4, show.plot = F, na.rm = F){

  lower_bound <- 10^(ndigits*-1)
  upper_bound <- 1 - lower_bound
  to_test <- seq(lower_bound, upper_bound, by = lower_bound)
  lik1 <- unlist(lapply (to_test, likelihood_pa, input, na.rm = na.rm))
  the_maxl <- to_test[which.max(lik1)]

  if (length (the_maxl) > 1){
    stop ('Poorly defined likelihood surface')
  }

  if (show.plot == T){
    plot (to_test,lik1,xlab = "Proportion of genes contributing to adaptation",ylab = "log Likelihood", type = "l")
    arrows (the_maxl,-10^50,the_maxl,10^50, col = "red", lty = 2)
  }

  the_maxl

}



