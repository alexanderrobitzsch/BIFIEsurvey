## File Name: BIFIE_lavaan_coef.R
## File Version: 0.05

BIFIE_lavaan_coef <- function(object)
{
    est <- object@Fit@x
    vcov1 <- BIFIE_lavaan_vcov(object=object)
    names(est) <- rownames(vcov1)
    # partable <- as.data.frame(object@ParTable)
    # res <- BIFIE_lavaan_lavInspect(object=object, what="est")
    return(est)
}
