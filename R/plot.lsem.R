## File Name: plot.lsem.R
## File Version: 0.367

#**** plot function for objects of class lsem
plot.lsem <- function( x, parindex=NULL, ask=TRUE, ci=TRUE,
            lintrend=TRUE, parsummary=TRUE,
            ylim=NULL, xlab=NULL, ylab=NULL, main=NULL, digits=3, ... )
{

    #--- change by Johanna Hartung 2016-06-20
    if (is.null(parindex)) {
        NP <- max(x$parameters$parindex)
        parindex <- 1:NP
    } else {
        NP <- parindex
    }
    NP <- max( parindex )

    ylim0 <- ylim
    if ( ! is.null(ylim) ){
        if ( ! is.list(ylim) ){
            ylim <- list()
            for (pp in 1:NP){
                ylim[[pp]] <- ylim0
            }
        }
    }
    xlab1 <- xlab
    if ( is.null(xlab)){
        xlab1 <- x$moderator
    }

    modgrid <- x$moderator.density

    #--- loop over all parameters
    for (pp in parindex ){
        graphics::par( mfrow=c(1,1))
        ind.pp <- which( parindex==pp)
        x.pp <- x$parameters
        x.pp <- x.pp[ x.pp$parindex==pp, ]
        t1 <- paste( x.pp$par[1] )
        if (! is.null(main) ){
            t1 <- main[ind.pp]
        }
        if (parsummary){
            t1 <- paste0( t1, "\n M=", round(x$parameters_summary$M[pp],digits),
                        ", SD=", round(x$parameters_summary$SD[pp],digits),
                        " | y="    , round(x$parameters_summary$lin_int[pp],digits)    ,
                        ifelse( x$parameters_summary$lin_slo[pp] < 0, "", "+" ),
                        round(x$parameters_summary$lin_slo[pp],digits),
                        "*m"
                                )
        }
        if ( is.null(ylim0) ){
            ylim1 <- NULL
            if (!is.null(x.pp$ci.lower)){
                ylim1 <- c( min( x.pp$ci.lower), max( x.pp$ci.upper) )
            }
        } else {
                ylim1 <- ylim[[ind.pp]]
        }
        ylim1_null <- is.null(ylim1)
        if (!ylim1_null){
            if (max(ylim1) %in% c("-Inf","Inf") ){
                ylim1_null <- TRUE
            }
        }
        if (ylim1_null){
            ylim1 <- range(x.pp$est)
        }
        moderator.grid <- x$moderator.grid
        xlim1 <- range( moderator.grid )
        if ( x$type=="MGM" ){
            moderator.grouped <- x$moderator.grouped
            xlim1 <- c( min( moderator.grouped$min), max(moderator.grouped$max ) )
        }
        ylab1 <- paste(x.pp$par[1])
        if ( ! is.null(ylab) ){
            ylab1 <- ylab[ind.pp]
        }

        x1 <- modgrid[,1]
        y1 <- x.pp$est

        graphics::plot( x1, y1, xlab=xlab1, ylab=ylab1,
                    main=t1, type="p", pch=16, ylim=ylim1, xlim=xlim1 )
        if ( x$type=="LSEM" ){
            graphics::lines( lsem_spline( x=x1, y=y1 )  )
        }
        if ( x$type=="MGM" ){
            G <- nrow(moderator.grouped)
            for (gg in 1:G){
                graphics::lines( moderator.grouped[gg,1:2], rep( x.pp$est[gg], 2 ) )
            }
            graphics::points( modgrid[,1], y=x.pp$est, pch=16 )
        }

        if (lintrend){
            x1a <- xlim1
            graphics::lines( x1a, x$parameters_summary$lin_int[pp] +
            x1a*x$parameters_summary$lin_slo[pp], lty=2, col=2)
        }
        if (ci & ( x.pp$op[1] !="fit" )  ){
            if ( x$type=="LSEM" ){
                graphics::lines( lsem_spline( x=x1, y=x.pp$ci.lower), lty=4)
                graphics::lines( lsem_spline( x=x1, y=x.pp$ci.upper), lty=4)
            }
            if ( x$type=="MGM" ){
                tick_length <- diff(xlim1) / 100
                for (gg in 1:G){
                    graphics::lines( rep( x1[gg], 2 ),
                                c( x.pp$ci.lower[gg], x.pp$ci.upper[gg]),
                                col="gray", lty=1)
                    graphics::points( x1[gg], x.pp$est[gg], pch=16 )
                    graphics::lines( x1[gg] + c(-1,1)*tick_length,  rep(x.pp$ci.lower[gg],2),
                                    col="gray", lty=1)
                    graphics::lines( x1[gg] + c(-1,1)*tick_length,  rep(x.pp$ci.upper[gg],2),
                                    col="gray", lty=1)
                }
            }
        }
        graphics::par(mfrow=c(1,1))
        graphics::par(ask=ask)
    }
    #***************** end loop
}

