## File Name: BIFIE_lavaan_vcov.R
## File Version: 0.05


BIFIE_lavaan_vcov <- function(object, ...)
{
    requireNamespace("lavaan")
    # res <- BIFIE_lavaan_lavInspect(object=object, what="vcov")
    # res <- as.matrix(res)
    lavaan_vcov <- methods::getMethod("vcov", "lavaan")
    vcov1 <- lavaan_vcov(object, ...)
    return(vcov1)
}
