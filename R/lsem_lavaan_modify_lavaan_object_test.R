## File Name: lsem_lavaan_modify_lavaan_object_test.R
## File Version: 0.02

lsem_lavaan_modify_lavaan_object_test <- function(object)
{
    # object@test
    test <- object@test
    test[[1]]$test <- "standard"
    object@test <- test

    #-- output
    return(object)
}
