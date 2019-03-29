## File Name: BIFIE_lavaan_summary.R
## File Version: 0.02

BIFIE_lavaan_summary <- function(object)
{
    requireNamespace("lavaan")
    lavaan_summary <- methods::getMethod("summary","lavaan")
    lavaan_summary(object)
}
