#' Calculate log-likelihood scores for estimate_pa.R
#'
#' @description \code{likelihood_pa} Calculates the log-likelihood scores
#' for use in estimate_pa.R. Setting na.rm = T assumes that missing values
#' are true negatives.
#'
#' @param input numeric vector
#' @param test_prop_adaptive numeric
#'
#' @return length-one numeric
#'
#' @examples
#' array1 <- array (0,c(5000,20))
#' array1[cbind(1:20,1:20)] <- 1
#' array1[1:5,1:3] <- 1
#' array1[6:10,4:6] <- 1
#' test_proportion <- 0.2
#' likelihood_pa (test_proportion, array1)
#'
#' @export


likelihood_pa <- function(test_prop_adaptive,input,na.rm = F){

  numcol <- ncol (input)
  adapted_per_gene <- rowSums (input,na.rm = na.rm)
  test_num_adaptive <- test_prop_adaptive * nrow(input)
  rate_adapted <- (sum(adapted_per_gene) / (test_num_adaptive * numcol))
  if (rate_adapted < 1){
    prob_per_gene <- (test_prop_adaptive * dbinom(adapted_per_gene,numcol,rate_adapted) ) + ((1 - test_prop_adaptive) * dbinom(adapted_per_gene,numcol, 0))
  } else {
    prob_per_gene <- NA
  }
  sum(log(prob_per_gene))

}



