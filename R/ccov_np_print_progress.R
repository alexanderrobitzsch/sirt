## File Name: ccov_np_print_progress.R
## File Version: 0.02
## File Last Change: 2020-02-14

ccov_np_print_progress <- function(progress, i, ii, display)
{
    if ( i < 20 ){
        if ( ii==display[i] & progress ){
            cat( paste( 5*i, "% ", sep="" ) )
            i <- i + 1
            if (i==11){
                cat("\n" )
            }
            utils::flush.console()
        }
    }
    return(i)
}
