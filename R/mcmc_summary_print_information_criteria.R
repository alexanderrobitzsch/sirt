## File Name: mcmc_summary_print_information_criteria.R
## File Version: 0.13
## File Last Change: 2018-12-30

mcmc_summary_print_information_criteria <- function(object)
{
    cat( "Number of iterations","=", object$iter, "\n" )
    cat( "Number of burnin iterations","=", object$burnin, "\n\n" )
    if ( ! is.null(object$description) ){
        cat( object$description, "\n")
    }
    cat("-----------------------------------------------------------------\n")
    cat( "Dbar","=", round( object$ic$Dbar, 2 ), "\n" )#, " | " )
    cat( "Dhat","=", round( object$ic$Dhat, 2 ), "\n" )#, " | " )
    cat( "pD","=", round( object$ic$pD, 2 ),  " | pD=Dbar - Dhat \n" )
    cat( "DIC","=", round( object$ic$DIC, 2 ),  " | DIC=Dhat + pD \n\n" )
}
