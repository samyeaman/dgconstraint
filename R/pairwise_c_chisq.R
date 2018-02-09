#' Calculate C-chisq
#'
#' @description \code{pairwise_c_chisq} returns the pairwise C-score using
#' the chi-square approach. Setting na.rm = T assumes that missing values
#' are true negatives.
#'
#' @param input numeric array
#' @param num_permute numeric
#' @param na.rm boolean
#'
#' @return length-one numeric
#'
#' @examples
#' array1 <- array (0,c(5000,20))
#' array1[cbind(1:20,1:20)] <- 1
#' array1[1:5,1:3] <- 1
#' array1[6:10,4:6] <- 1
#' pairwise_c_chisq (array1)
#'
#' @export

pairwise_c_chisq <- function (input,num_permute = 10000,na.rm = F){

  numcol <- ncol (input)
  c_score <- array(NA, c(numcol,numcol))

  for (loop1 in 1:(numcol-1)){
    for (loop2 in (loop1+1):numcol){
      results_chisq <- array (NA,num_permute)
      input_sub <- as.matrix (input[,c(loop1,loop2)])
      obs1 <- rowSums (input_sub,na.rm = na.rm)
      exp1 <- array (mean(obs1), length(obs1))
      chisq1 <- (obs1 - exp1)^2 / exp1

      for (loop3 in 1:num_permute){
        input_sub[,1] <- sample (input_sub[,1],nrow (input_sub),replace = F)
        input_sub[,2] <- sample (input_sub[,2],nrow (input_sub),replace = F)

        obs2 <- rowSums(input_sub,na.rm = na.rm)
        exp2 <- array (mean(obs2),length(obs2))
        chisq2 <- (obs2 - exp2)^2/exp2

        results_chisq[loop3] <- sum (chisq2)
      }
      c_score [loop1,loop2]  <- (sum (chisq1) - mean (results_chisq)) / sd (results_chisq)
    }
  }

  mean(c_score,na.rm = T)

}
