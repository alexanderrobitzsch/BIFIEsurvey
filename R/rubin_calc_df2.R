## File Name: rubin_calc_df2.R
## File Version: 0.05



######################################################
rubin_calc_df2 <- function( B, W, Nimp, indices=NULL, digits=2)
{
    if ( ! is.null(indices) ){
        W <- W[indices]
        B <- B[indices]
    }
    B <- B + W * 1E-15
    df <- ( 1 + Nimp / ( Nimp + 1 ) * W / B  )^2 * (Nimp - 1 )
    df <- round( df, digits )
    df <- ifelse( ( df > 1000 ) | ( Nimp==1 ), Inf, df )
    return(df)
}
######################################################
