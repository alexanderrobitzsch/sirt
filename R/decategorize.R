## File Name: decategorize.R
## File Version: 0.122
## File Last Change: 2023-03-10

#* decategorize
decategorize <- function( dat, categ_design=NULL )
{
    # preliminaries
    dat4 <- dat3 <- dat
    dfr <- categ_design

    #** handle categories
    if ( ! is.null( dfr ) ){
        vars <- sort( unique( paste( dfr$variable )))
        VV <- length(vars)
        for (vv in 1:VV){
            dfr.vv <- dfr[ paste(dfr$variable)==vars[vv], ]
            dat4[, vars[vv] ] <- dfr.vv[ match( dat3[,vars[vv]], dfr.vv$recode ), 'orig']
        }
    }
    return(dat4)
}

