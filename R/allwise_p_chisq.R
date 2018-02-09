#' Calculate p-value for C-chisq
#'
#' @description \code{allwise_p_chisq} returns the p-value for 'all lineages'
#' contrast using chi-square. Setting na.rm = T assumes that missing values
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
#' allwise_p_chisq (array1)
#'
#' @export

allwise_p_chisq <- function (input,num_permute = 10000,na.rm = F){

  results_chisq <- array (NA,num_permute)

  obs1 <- rowSums (input,na.rm = na.rm)
  exp1 <- array (mean(obs1), length(obs1))
  chisq1 <- (obs1 - exp1)^2 / exp1

  for (loop1 in 1:num_permute){

    input_sub <- input

    for (loop2 in 1:ncol(input_sub)){
      input_sub[,loop2] <- sample (input_sub[,loop2],nrow (input_sub),replace = F)
    }

    obs2 <- rowSums(input_sub,na.rm = na.rm)
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
