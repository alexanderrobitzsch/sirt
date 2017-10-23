## File Name: rm_facets_print_progress_trait_distribution.R
## File Version: 0.01
## File Last Change: 2017-10-02 15:42:29


rm_facets_print_progress_trait_distribution <- function( parm, parmlabel , digits_trait )
{
    cat( paste( " ", parmlabel  , "= " , round( parm , digits_trait ) , sep="") , "\n")	
}
