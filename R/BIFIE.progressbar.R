## File Name: BIFIE.progressbar.R
## File Version: 0.09


#--- Computation of a progress bar
BIFIE.progressbar <- function( ops, prblen )
{
    prb <- prblen
    vec <- seq(1, ops)
    vec[ ops ] <- ops - 0.1
    NR <- ops / prb
    m1 <- vec %% NR
    pr1 <- 1 * ( diff(m1) < 0  )
    pr1 <- c( 1, pr1 )
    # returns a vector of zeroes and one indicating
    # iteration of a move in th progress bar
    return(pr1)
}
