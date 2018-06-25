## File Name: wle.rasch.jackknife.R
## File Version: 1.07


###############################################################################################
# Jackknife standard error estimation of WLE
wle.rasch.jackknife <- function( dat, b, itemweights=1+0*b, pid=NULL,
            testlet=NULL, stratum=NULL,
            size.itempop=NULL )
{
        # INPUT:
        # size.itempop ... finite sampling correction
        #
        #*************************************************
        # define person ID
        if( is.null(pid) ){ pid <- seq( 1, nrow(dat)) }
        # WLE estimate
        wle <- wle.rasch( dat=dat, b=b, itemweights=itemweights )
        dat.resp <- 1-is.na( dat )
        dfr <- est.bias <- NULL

        if ( ! is.null( size.itempop ) ){
            if ( is.null( names( size.itempop ) ) ){
                    names(size.itempop) <- sort( unique( stratum ) )
                            }}

        cat("Jackknife standard of WLE\n")
        #*********************************
        # stratified item sampling
        if ( ! is.null(stratum)){
            strata <- unique(stratum)
            cat("Stratified item sampling\n")
            S <- length(strata)
            wlevarstrata <- matrix( 0, nrow=nrow(dat), ncol=S )
            colnames(wlevarstrata) <- strata
            items.in.strata <- wlevarstrata
            cat(paste(S, "strata\n") )
            # save data response matrix
            dat.resp <- wle$dat.resp
            sk <- 0
            for ( ss in strata ){
#                ss <- strata[3]
                items.ss <- colnames(dat)[ stratum==ss ]
                II <- length(items.ss)
                wlejack <- matrix( 0, nrow=nrow(dat), ncol=II)
                for (ii in 1:II){   # begin loop: items within a stratum
#                    ii <- 1
                    # calculate item index
                    ind.ii <- which( colnames(dat)==items.ss[ii] )
                    # calculate weight matrix of stratification
                    ind.ss <- which( stratum==ss )
                    dat.resp.ss <- dat.resp[, -ind.ii ]
                    # sum of itemweights within a stratum
                    sumitemweights <- rowSums( dat.resp[, ind.ss ] )
                    ind1 <- intersect( colnames(dat.resp.ss), items.ss )
                    h1 <- as.matrix( dat.resp.ss[, ind1 ], ncol=length(ind1))
                    dat.resp.ss[, ind1] <- sumitemweights / rowSums( h1 ) * dat.resp.ss[, ind1 ]
                    dat.resp.ss[ is.na(dat.resp.ss) ] <- 0
                    # WLE estimation
                    wlejack[,ii] <- wle.rasch( dat=dat[, - ind.ii ], dat.resp=dat.resp.ss,
                                            b=b[ - ind.ii  ],
                                            theta=wle$theta  )$theta
                                }       # end loop: items within a stratum
                    ind.resp.ss <- 1*( dat.resp[, items.ss ] > 0 )
                    wlejack <- wlejack * ind.resp.ss
                    N.wlejack <- rowSums( ind.resp.ss )
                    items.in.strata[,ss] <- N.wlejack

                    mean.wlejack <- rowSums( wlejack ) / N.wlejack
                    i1 <- rowSums( ( wlejack - outer( mean.wlejack, rep(1, ncol(wlejack) ) ) )^2 )
                    # finite sampling correction
                    if ( is.null(size.itempop[ss]) ){ fcs <- rep(1,nrow(dat)) } else
                                {
                                    fcs <- 1 - N.wlejack / size.itempop[ss]
                            cat( "\nStratum ", ss, "(Mean) Correction Factor", round(mean(fcs),5), "\n")
                                }
                    wlevarstrata[,ss] <- fcs * (N.wlejack -1 )/ N.wlejack  * i1
                    sk <- sk+1
                    cat(paste(sk, ".",sep="")) ;
                    utils::flush.console() ; if( sk %% 10==0 | sk==S ){ cat("\n") }
                        }   # end loop strata
            jack.se <- sqrt( rowSums(wlevarstrata )   )
            jackunits <- NA
                    }
        #*********************************
        # random sampling of testlets
        if ( ! is.null(testlet) ){
            testlets <- unique(testlet)
            cat("Simple random sampling of testlets\n")
            I <- length(testlets)
            wlejack <- matrix( 0, nrow=nrow(dat), ncol=I )
            testlet.ind <- wlejack
            cat(paste(I, "testlets\n") )
            testletcount <- rep(0,nrow(dat))
                for (ii in 1:I){
                    # ii <- 1
                    wlejack[,ii] <- wle.rasch( dat=dat[, testlet !=testlets[ii] ], b=b[ testlet !=testlets[ii]  ],
                                            itemweights=itemweights[ testlet !=testlets[ii]  ],
                                            theta=wle$theta  )$theta
                    testletcount <- testletcount + ( rowSums( wle$dat.resp[, testlet==testlets[ii] ] ) > 0 )
                    testlet.ind[,ii] <- ( rowSums( wle$dat.resp[, testlet==testlets[ii] ] ) > 0 )
                    cat(paste(ii, ".",sep="")) ;
                    utils::flush.console() ; if( ii %% 10==0 | ii==I ){ cat("\n") }
                            }
                # number of Jackknife units
                jackunits <- testletcount
            # finite sampling correction
            N.items <- rowSums( wle$dat.resp > 0 )
            if ( is.null(size.itempop ) ){ fcs <- rep(1,nrow(dat)) } else
                    { fcs <- 1 - N.items / size.itempop }
        # include indicator matrix of testlet responses
                jack.se <- sqrt( fcs * (jackunits - 1) / jackunits * rowSums( testlet.ind*( wlejack - outer( wle$theta, rep(1,I) ) )^2  ) )
                        }
        #*********************************
        # Simple random item sampling
        if ( is.null(testlet) & ( is.null(stratum) ) ){     # begin simple random sampling
            I <- ncol(dat)
            wlejack <- matrix( 0, nrow=nrow(dat), ncol=I )
            cat("Simple random sampling of items\n")
            cat(paste(I, "items\n") )
            for (ii in 1:I){
                # ii <- 1
                wlejack[,ii] <- wle.rasch( dat=dat[, -ii ], b=b[-ii], itemweights=itemweights[-ii],
                            theta=wle$theta  )$theta
                cat(paste(ii, ".",sep="")) ;
                utils::flush.console() ; if( ii %% 10==0 | ii==I ){ cat("\n") }
                        }
            # estimate WLE bias
            dfr <- data.frame( "wle"=wle$theta )
            UJ <- rowSums( dat.resp )
            mean.jack <- rowSums( wlejack, na.rm=T ) / UJ
            dfr$wle.jack <- UJ * wle$theta - (UJ-1)*mean.jack
            est.bias <- dfr$est.bias <- - dfr$wle.jack + dfr$wle
            # number of Jackknife units
            jackunits <- rowSums( wle$dat.resp > 0 )
            # finite sampling correction
            if ( is.null(size.itempop ) ){ fcs <- rep(1,nrow(dat)) } else
                    { fcs <- 1 - jackunits / size.itempop }
            jack.se <- sqrt( fcs * (jackunits - 1) / jackunits * rowSums( ( wle$dat.resp>0)*( wlejack - outer( wle$theta, rep(1,I) ) )^2  ) )
                            }   # end simple random sampling
        #*********************************
        # calculate and print reliability
        v1 <- stats::var(wle$theta)
        v2 <- mean( jack.se^2    )
        wle.rel <-     ( v1 - v2 ) / v1
        cat("WLE Reliability=", round(wle.rel,3), "\n")
        #*********************************
        # collect results
        dfr <- data.frame( "pid"=pid, "N.jackunits"=jackunits,
                            "sum.itemweights"=rowSums( wle$dat.resp),
                            "score"=rowSums( dat, na.rm=T), "max"=rowSums( 1 -is.na(dat) ),
                            "wle"=wle$theta, "wle.jackse"=jack.se,
                            "wlejack"=wlejack,
                            "est.bias"=ifelse( is.null( est.bias), NA, est.bias )
                            )
        res <- list( "wle"=dfr, "wle.rel"=wle.rel )
        return(res)
}
###############################################################################################

