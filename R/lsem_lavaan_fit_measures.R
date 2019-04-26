## File Name: lsem_lavaan_fit_measures.R
## File Version: 0.07

lsem_lavaan_fit_measures <- function(object, fit_measures)
{
    fM <- sirt_import_lavaan_fitMeasures(object=object, fit_measures=fit_measures)
    fit_measures <- intersect( fit_measures, names(fM))
    fM <- fM[ fit_measures ]
    return(fM)
}
