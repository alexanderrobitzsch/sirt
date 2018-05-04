## File Name: sirt_summary_print_packages.R
## File Version: 0.02

sirt_summary_print_packages <- function(packages)
{
    for (pack in packages){
        sirt_summary_print_package(pack=pack)
    }
}
