## File Name: package_version_date.R
## File Version: 0.08
## File Last Change: 2018-12-30

package_version_date <- function(package)
{
    d1 <- utils::packageDescription(pkg=package)
    res <- paste( d1$Package, " ", d1$Version,
                " (", d1$Date, ")", sep="")
    return(res)
}
