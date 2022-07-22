#' Calculate p-value for C-hyper
#'
#' @description \code{single_p_hyper} returns the p-value for a single pair
#' of lineages using the hypergeometric approach. Setting na.rm = T assumes 
#' that missing values are true negatives.
#'
#' @param input1 numeric vector
#' @param input2 numeric vector
#' @param na.rm boolean
#'
#' @return length-one numeric
#'
#' @examples
#' vector1 <- rep (0,1000)
#' vector2 <- rep (0,1000)
#' vector1[c(1:5,15:20,24,31:35,47:50,75)] <- 1
#' vector2[c(6:9,12,19:23,28:31,78:80)] <- 1
#' single_p_hyper (vector1, vector2)
#'
#' @export


single_p_hyper <- function (input1, input2, na.rm = F){
  
  ax <- sum (input1, na.rm = na.rm)
  ay <- sum (input2, na.rm = na.rm)
  g0 <- as.numeric (length (input1))

  sd_hyp <- sqrt((ax*ay)*(g0-ax)*(g0-ay)/(g0^2*(g0-1)))
  exp_hyp <- ax * ay / g0
  obs_hyp <- sum (input1 == 1 & input2 == 1, na.rm = na.rm)

        if (sd_hyp != 0){
          result <- (obs_hyp - exp_hyp) / sd_hyp
        } else {
          result <- 0
          warning ('Pair of lineages has no shared adapted loci')
        }
  
  min <- min(ax,ay)
  max <- max(ax,ay) 
  p <- sum(dhyper(obs_hyp:min, min, g0 - min, max))
  p

}
