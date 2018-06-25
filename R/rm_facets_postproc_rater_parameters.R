## File Name: rm_facets_postproc_rater_parameters.R
## File Version: 0.07

rm_facets_postproc_rater_parameters <- function( rater.index, dat2, dat2.resp, b.rater, a.rater,
        rater.index1, rater_item_int )
{
    N <- colSums( dat2.resp )
    M1 <- colSums( dat2 ) / N
    N <- stats::aggregate( N, list( rater.index ), sum, na.rm=TRUE )[,2]
    M1 <- stats::aggregate( M1, list( rater.index ), mean, na.rm=TRUE )[,2]
    rater <- data.frame( "rater"=rater.index1[,1], "N"=N, "M"=M1,     "b"=b.rater,    "a"=a.rater )
    rater$thresh <- rater$a * rater$b
    if ( ! rater_item_int){
        rater$b.cent <- rm_facets_center_value(x=rater$b, value=0)
        rater$a.cent <- rm_facets_center_value(x=rater$a, value=1)
    }
    if ( rater_item_int){
        vec <- rm_facets_string_part_extract( x=rater$rater, split="-", part=2)
        rater$b.cent <- rm_facets_center_value_aggregate(x=rater$b, index=vec, value=0)
        rater$a.cent <- rm_facets_center_value_aggregate(x=rater$a, index=vec, value=1)
    }
    #--- output
    return(rater)
}

