## File Name: package_version_date.R
## File Version: 0.02
## File Last Change: 2017-02-17 13:02:52

package_version_date <- function(package)
{
    d1 <- utils::packageDescription(pkg=package)
	res <- paste( d1$Package , " " , d1$Version ,
				" (" , d1$Date , ")" , sep="")	
	return(res)
}
