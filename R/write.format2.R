## File Name: write.format2.R
## File Version: 1.09



#----------------------------------------------------------------
# utility function for formatting output in write.fwf2
write.format2 <- function( vec1, ff, fr ){
    if (fr==0){
        vec2 <- round( vec1, fr )
        blank.vv <- paste( rep( " ", ff ), collapse="" )
        vec2 <- paste( substring( blank.vv, 1, ff - nchar(vec2) ), vec2, sep="")
            } else {
        d.vv <- round( vec1, fr ) + 10^(-(fr+1))
        # generate blank
        blank.vv <- paste( rep( " ", ff+1 ), collapse="" )
        d.vv <- paste( substring( blank.vv, 1, ff+1 - nchar(d.vv) ), d.vv, sep="")
        g.vv <- grep("NA",d.vv)
        d.vv[ g.vv  ] <- ifelse( ff > 1,  gsub( "NA", " .", d.vv[g.vv] ), gsub( "NA", ".", d.vv[g.vv] ) )
        vec2 <- substring( d.vv, 1, ff )
        vec2
            }
    return(vec2 )
    }
#---------------------------------------------------------------

.write.format2 <- write.format2
