## File Name: rm_facets_center_value.R
## File Version: 0.01

rm_facets_center_value <- function(x, value=0)
{
	y <- x + value - mean(x, na.rm=TRUE)
	return(y)
}
	
