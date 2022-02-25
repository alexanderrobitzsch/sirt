## File Name: mgsem_proc_model_include_missing_entries.R
## File Version: 0.09


mgsem_proc_model_include_missing_entries <- function(model_hh, types,
        entries=c("est","index"), I, D)
{
    for (tt in types){
        val <- 0
        if (tt=="NU"){ NR <- I; NC <- 1 }
        if (tt=="ALPHA"){ NR <- D; NC <- 1 }
        if (tt=="LAM"){ NR <- I; NC <- D }
        if (tt=="PHI"){ NR <- D; NC <- D }
        if (tt=="PSI"){ NR <- I; NC <- I }
        if (tt=="B"){ NR <- D; NC <- D }
        for (ee in entries){
            if (is.null(model_hh[[ee]][[tt]])){
                model_hh[[ee]][[tt]] <- matrix(val, nrow=NR, ncol=NC)
            }
        }
    }
    return(model_hh)
}
