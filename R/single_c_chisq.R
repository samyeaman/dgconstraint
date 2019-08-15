#' Calculate C-chisq
#'
#' @description \code{single_c_chisq} returns the C-score for a given 
#' pair of lineages, using the chi-square approach. Setting na.rm = T 
#' assumes that missing values are true negatives.
#'
#' @param input1 numeric vector
#' @param input2 numeric vector
#' @param num_permute numeric
#' @param na.rm boolean
#'
#' @return length-one numeric
#'
#' @examples
#' vector1 <- rep (0,1000)
#' vector2 <- rep (0,1000)
#' vector1[c(1:5,15:20,24,31:35,47:50,75)] <- 1
#' vector2[c(6:9,12,19:23,28:31,78:80)] <- 1
#' single_c_chisq (vector1, vector2)
#'
#' @export

single_c_chisq <- function (input1, input2 ,num_permute = 10000,na.rm = F){

  results_chisq <- array (NA,num_permute)
  
  if (length(input1) != length(input2)){
  	warning('Lineages are not the same length')
  }	
  
  both <- cbind(input1, input2)
  obs1 <- rowSums (both,na.rm = na.rm)
  exp1 <- array (mean(obs1), length(obs1))
  chisq1 <- (obs1 - exp1)^2 / exp1

  for (loop3 in 1:num_permute){
    sub1 <- sample(input1,length (input1),replace = F)
    sub2 <- sample(input2, length(input2), replace = F)
	both2 <- cbind(sub1, sub2)
    obs2 <- rowSums(both2,na.rm = na.rm)
    exp2 <- array (mean(obs2),length(obs2))
    chisq2 <- (obs2 - exp2)^2/exp2

    results_chisq[loop3] <- sum (chisq2)
  }
      
  c_score <- (sum (chisq1) - mean (results_chisq)) / sd (results_chisq)
  c_score
}
