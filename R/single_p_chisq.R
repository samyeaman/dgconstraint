#' Calculate p-value for C-chisq
#'
#' @description \code{single_p_chisq} returns the p-value for a single pair
#' of lineages, contrasted using chi-square. Setting na.rm = T assumes that 
#' missing values are true negatives.
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
#' single_p_chisq (vector1, vector2)
#'
#' @export

single_p_chisq <- function (input1, input2 ,num_permute = 10000,na.rm = F){

  results_chisq <- array (NA,num_permute)
  
  if (length(input1) != length(input2)){
  	warning('Lineages are not the same length')
  }	

  both <- cbind(input1, input2)
  obs1 <- rowSums (both,na.rm = na.rm)
  exp1 <- array (mean(obs1), length(obs1))
  chisq1 <- (obs1 - exp1)^2 / exp1

  for (loop1 in 1:num_permute){

    in1 <- sample(input1, length(input1), replace = F)
    in2 <- sample(input2, length(input2), replace = F)
    
    both2 <- cbind(in1, in2)
    obs2 <- rowSums(both2,na.rm = na.rm)
    exp2 <- array (mean(obs2),length(obs2))
    chisq2 <- (obs2 - exp2)^2/exp2

    results_chisq[loop1] <- sum (chisq2)

  }

  out1 <- sum (results_chisq > sum(chisq1))

  if (out1 == 0){
    warning('Note: estimated p-value is lower than 1/(number of permutations)')
    1/num_permute
  } else {
    out1/num_permute
  }

}
