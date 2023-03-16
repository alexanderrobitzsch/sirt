## File Name: sirt_summary_print_package.R
## File Version: 0.06
## File Last Change: 2023-03-08

sirt_summary_print_package <- function(pack)
{
    cat( package_version_date(package=pack), '\n' )
}
