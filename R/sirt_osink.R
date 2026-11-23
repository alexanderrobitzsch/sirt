## File Name: sirt_osink.R
## File Version: 0.09


sirt_osink <- function(file)
{
    CDM::osink( file=file, suffix=paste0( '__SUMMARY.Rout') )
}
