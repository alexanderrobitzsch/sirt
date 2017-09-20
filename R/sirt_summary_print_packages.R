## File Name: sirt_summary_print_packages.R
## File Version: 0.02
## File Last Change: 2017-09-20 10:51:21

sirt_summary_print_packages <- function(packages)
{
	for (pack in packages){
		sirt_summary_print_package(pack=pack)
	}
}
