## File Name: xxirt_print_progress.R
## File Version: 0.072

xxirt_print_progress <- function(dev, dev0, dev00, pen_val, conv0, conv1,
        iter, verbose1, verbose2, verbose_index, verbose3=FALSE,
        opt_fun_value, prior_par )
{
    if (verbose1){
        cat( paste( '   Optimization function value=', round( opt_fun_value, 4 ),
            if (iter > 1 ){ ' | Change=' } else {''},
            if( iter > 1){ round( - dev + dev0, 6 )} else { ''}    ,'\n',sep='') )
        cat( paste( '   -Log-likelihood=',
                paste( round( dev00/2, 4), collapse=' ' ), '\n', sep=''))
        cat( paste( '   Prior function value=',
                paste( round( prior_par, 6), collapse=' ' ), '\n', sep=''))
        cat( paste( '   Penalty function=',
                paste( round( pen_val, 6), collapse=' ' ), '\n', sep=''))
        cat( paste( '   Deviance=',
                paste( round( dev00, 4), collapse=' ' ), '\n', sep=''))
        cat( paste( '    Maximum item parameter change=',
                paste( round( conv0, 6), collapse=' ' ), '\n', sep=''))
        cat( paste( '    Maximum theta distribution parameter change=',
                paste( round( conv1, 6), collapse=' ' ), '\n', sep=''))
        utils::flush.console()
    }
    if (verbose2){
        v1 <- paste0( verbose_index, ' | ',
                'Iteration=', iter, ' | Optimization function value=',
                round(dev,4), '\n')
        cat(v1)
        utils::flush.console()
    }
    if (verbose3){
        cat('-')
        utils::flush.console()
    }
}
