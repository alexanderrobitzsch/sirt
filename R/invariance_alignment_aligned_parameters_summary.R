## File Name: invariance_alignment_aligned_parameters_summary.R
## File Version: 0.07

invariance_alignment_aligned_parameters_summary <- function(x, label=NULL)
{
    dfr <- data.frame(Med=sirt_colMedians(x=x), M=colMeans(x=x, na.rm=TRUE),
                        SD=sirt_colSDs(x=x), Min=sirt_colMins(x=x),
                        Max=sirt_colMaxs(x=x) )
    if ( ! is.null(label) ){
        colnames(dfr) <- paste0( colnames(dfr), '.', label )
    }
    return(dfr)
}
