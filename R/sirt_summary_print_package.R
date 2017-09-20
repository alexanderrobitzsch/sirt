## File Name: sirt_summary_print_package.R
## File Version: 0.02
## File Last Change: 2017-09-20 10:49:09

sirt_summary_print_package <- function(pack)
{
	cat( package_version_date(package=pack) , "\n" )
}
