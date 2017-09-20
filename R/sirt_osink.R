## File Name: sirt_osink.R
## File Version: 0.01
## File Last Change: 2017-09-20 10:37:43


sirt_osink <- function(file)
{
    CDM::osink( file = file, suffix = paste0(  "__SUMMARY.Rout") )
}
