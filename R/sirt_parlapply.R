## File Name: sirt_parlapply.R
## File Version: 0.02

sirt_parlapply <- function(cl, X, FUN, verbose=FALSE, ...)
{
    args <- list(...)
    args$cl <- cl
    args$X <- X
    if (verbose){
        what <- pbapply::pblapply
        args$FUN <- FUN
    } else {
        what <- parallel::parLapply
        args$fun <- FUN
    }
    res_all <- do.call(what=what, args=args)
    return(res_all)
}
