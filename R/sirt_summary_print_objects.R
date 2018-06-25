## File Name: sirt_summary_print_objects.R
## File Version: 0.09

sirt_summary_print_objects <- function(obji, from=NULL, to=NULL,
        digits=3, rownames_null=TRUE)
{
    #------ data frame or matrix
    if ( is.data.frame(obji) | is.matrix(obji)){
        if (is.null(from)){
            from <- 1
        }
        if (is.null(to)){
            to <- ncol(obji)
        }
        if ( length(digits)==1){
            rvars <- seq( from, to)
            digits <- sirt_vector_with_names(value=digits, names=colnames(obji)[rvars] )
        }
        if ( length(digits) > 1){
            rvars <- names(digits)
        }
        for (vv in rvars ){
            obji[,vv] <- sirt_round_vector(x=obji[,vv], digits=digits[vv])
        }
        if (rownames_null){
            rownames(obji) <- NULL
        }
    }
    #------ vector
    if ( is.vector(obji) ){
        obji <- sirt_round_vector( x=obji, digits=digits )
    }
    #--- print object
    print(obji)
}
