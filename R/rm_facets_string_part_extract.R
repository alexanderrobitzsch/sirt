## File Name: rm_facets_string_part_extract.R
## File Version: 0.06

rm_facets_string_part_extract <- function( x, split, part)
{
    vec <- strsplit( paste(x), split=split )
    vec <- unlist( lapply( vec, FUN=function(vv){ vv[part] } ) )
    return(vec)
}
