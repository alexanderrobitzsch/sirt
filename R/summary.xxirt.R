## File Name: summary.xxirt.R
## File Version: 0.275


#--- summary for xxirt object
summary.xxirt <- function( object, digits=3, file=NULL, ...)
{
    # open sink
    sirt_osink( file=file )

    # print summary
    res <- xxirt_summary_parts(object=object, digits=digits)

    # close sink
    sirt_csink( file=file )
}
