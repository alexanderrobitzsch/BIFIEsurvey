## File Name: rubin_calc_df.R
## File Version: 0.10

######################################################
# calculate degrees of freedom according to Rubin
rubin_calc_df <- function( res_pars, Nimp, indices=NULL, digits=2)
{
    W <- res_pars$pars_varWithin
    B <- res_pars$pars_varBetween
    df <- rubin_calc_df2( B=B, W=W, Nimp=Nimp, indices=indices, digits=digits)
    return(df)
}
######################################################
