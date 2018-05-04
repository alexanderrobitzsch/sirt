## File Name: summary.lsdm.R
## File Version: 0.01


#.............................
# summary for LSDM function   
summary.lsdm <- function( object , ... )
{
        lsdmobj <- object
        # generate sequence for display
        display.separate <- paste( rep("." , each=80 ) , collapse="" )
        # display progress
        cat( display.separate , "\n" )
        cat( "LSDM -- Least Squares Distance Method of Cognitive Validation \n")
        cat("Reference: Dimitrov, D. (2007) Applied Psychological Measurement, 31, 367-387.\n")
        cat( display.separate , "\n" ) ; flush.console()
        cat("\nModel Fit\n\n")
        cat(paste( "Model Fit LSDM   -  Mean MAD:" , formatC( round( lsdmobj$mean.mad.lsdm0 , 3 ),digits=3 , width=6) , 
                        "    Median MAD:" , formatC( round( median(lsdmobj$item$mad.lsdm) , 3 ),digits=3 , width=6)     , "\n") )
        cat(paste( "Model Fit LLTM   -  Mean MAD:" , formatC( round( lsdmobj$mean.mad.lltm , 3 ),digits=3, width=6) , 
                    "    Median MAD:" , formatC( round( median(lsdmobj$item$mad.lltm) , 3 ),digits=3 , width=6) ,        
                    "   R^2=" , format( round( summary(lsdmobj$lltm)$r.squared , 3 ),digits=3) ,   "\n") )
        cat( display.separate , "\n" )
        cat("\nAttribute Parameters\n\n")
        dfr.a <- data.frame( "N.Items" = colSums(lsdmobj$Qmatrix) , round( lsdmobj$attr.pars , 3 )  )
        elim.a <- union( grep("Q" , colnames(dfr.a) ) , grep( "sigma" , colnames(dfr.a) ) )
        print( dfr.a[ , -elim.a ] )
        cat( display.separate , "\n" )
        cat("\nItem Parameters\n\n")
        dfr.i <- round( lsdmobj$item , 3 )

        displ.i <- sort( union(  which( colnames(dfr.i)   %in% c(  "a.2PL" , "b.2PL" , "N.Items" , "b.1PL" ) )  ,
                            grep( "mad" , colnames(dfr.i)   ) )                        
                                        )
        print( round( dfr.i[ ,  displ.i  ] , 3 ) )
        cat( display.separate , "\n" )
        cat("\nDiscrimination Parameters\n\n")
        dfr.i <- round( lsdmobj$W , 3 )
        print( dfr.i )

}


