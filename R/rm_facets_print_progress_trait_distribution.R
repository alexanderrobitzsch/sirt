## File Name: rm_facets_print_progress_trait_distribution.R
## File Version: 0.06


rm_facets_print_progress_trait_distribution <- function( parm, parmlabel, digits_trait )
{
    cat( paste( " ", parmlabel, "=", round( parm, digits_trait ), sep=""), "\n")
}
