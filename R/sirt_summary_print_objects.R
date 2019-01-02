## File Name: sirt_summary_print_objects.R
## File Version: 0.19

sirt_summary_print_objects <- function(obji, from=NULL, to=NULL,
        digits=3, rownames_null=TRUE, grep_string=NULL)
{
    #------ data frame or matrix
    if ( is.data.frame(obji) | is.matrix(obji)){
        include_names <- FALSE
        if ( is.null(colnames(obji) ) ){
            include_names <- TRUE
            colnames(obji) <- paste0("V",seq_len(ncol(obji)) )
        }
        if ( ! is.null(grep_string) ){
            obji <- obji[, grep( grep_string, colnames(obji)) ]
        }
        if (is.null(from)){
            from <- 1
        }
        if (is.null(to)){
            to <- ncol(obji)
        }
        if ( length(digits)==1){
            rvars <- seq( from, to)
            names1 <- colnames(obji)[rvars]
            if (is.null(names1)){
                names1 <- rvars
            }
            digits <- sirt_vector_with_names(value=digits, names=names1 )
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
        if (include_names){
            colnames(obji) <- NULL
        }
    }
    #------ vector
    if ( is.vector(obji) ){
        obji <- sirt_round_vector( x=obji, digits=digits )
    }
    #--- print object
    print(obji)
}
