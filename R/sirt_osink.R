## File Name: sirt_osink.R
## File Version: 0.07


sirt_osink <- function(file)
{
    CDM::osink( file=file, suffix=paste0(  '__SUMMARY.Rout') )
}
