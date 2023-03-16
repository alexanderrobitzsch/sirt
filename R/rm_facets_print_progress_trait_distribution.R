## File Name: rm_facets_print_progress_trait_distribution.R
## File Version: 0.08
## File Last Change: 2018-12-30


rm_facets_print_progress_trait_distribution <- function( parm, parmlabel, digits_trait )
{
    cat( paste( " ", parmlabel, "=", round( parm, digits_trait ), sep=""), "\n")
}
