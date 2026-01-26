## File Name: linking_haebara_proc_estimated_dif_effects.R
## File Version: 0.02

linking_haebara_proc_estimated_dif_effects <- function(mat)
{
    ind <- which(rowSums(1-is.na(mat))==1)
    if (length(ind)>0){
        mat[ind,] <- NA
    }
    return(mat)
}
