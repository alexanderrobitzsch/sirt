## File Name: sirt_osink.R
## File Version: 0.01


sirt_osink <- function(file)
{
    CDM::osink( file = file, suffix = paste0(  "__SUMMARY.Rout") )
}
