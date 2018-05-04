## File Name: lsem_fitsem.R
## File Version: 0.36

##############################################################
lsem_fitsem <- function( dat , weights , lavfit ,
            fit_measures , NF , G , moderator.grid , verbose ,
            pars , standardized, variables_model,
            sufficient_statistics , lavaan_fct , lavmodel , 
            ... )
{

    parameters <- NULL
    fits <- NULL     
    pars0 <- pars    
    env_temp <- environment()
    
    if (verbose){
        cat( "** Fit lavaan model\n")
        G1 <- min(G,10)    
        pr <- round( seq(1,G , len=G1) )
        cat("|")
        cat( paste0( rep("*",G1) , collapse="") )
        cat("|\n")
        cat("|")
    }

    for (gg in 1:G){
        # gg <- 1
        dat$weight <- weights[,gg]
        #***** fit the model using weighted data
        if (! sufficient_statistics){
            datsvy <- survey::svydesign(id=~index, weights=~weight, data=dat)        
            # assign(x="lavmodel__", value=lavmodel, pos=1)
            assign_args <- list( x="lavmodel__", value=lavmodel, pos=1)
            res0 <- do.call( what="assign", args=assign_args)
            survey.fit <- lavaan.survey::lavaan.survey(lavaan.fit=lavfit, 
                                survey.design=datsvy )
        }                        
        #***** fit the model using sufficient statistics
        if (sufficient_statistics){
            res <- lsem_weighted_mean( x=dat[ , variables_model], weights=dat$weight )
            wmean <- res$mean
            res <- lsem_weighted_cov( x=dat[ , variables_model], weights=dat$weight )
            wcov <- res$cov
            Nobs <- round( res$Nobs )
            if (lavaan_fct=="sem"){
                survey.fit <- lavaan::sem(model = lavmodel, sample.cov = wcov , 
                                sample.mean = wmean , sample.nobs = Nobs , ... )
            }
            if (lavaan_fct=="lavaan"){
                survey.fit <- lavaan::lavaan(model = lavmodel, sample.cov = wcov , 
                                sample.mean = wmean , sample.nobs = Nobs , ... )
            }
        }    
                            
        dfr.gg <- pars <- lavaan::parameterEstimates(survey.fit)                 
        if (standardized){            
            sol <- lavaan::standardizedSolution( survey.fit )
            colnames(sol)[ which( colnames(sol) == "est.std" ) ] <- "est"
            sol$lhs <- paste0( "std__" , sol$lhs)
            pars <- sirt_rbind_fill( x=pars, y=sol )
            # pars <- plyr::rbind.fill( pars , sol )    
            dfr.gg <- pars
        }                     
        pars <- paste0( pars$lhs , pars$op , pars$rhs )                    
        NP <- length(pars0)
        ind <- match( pars0 , pars )
        dfr.gg <- dfr.gg[ ind , ]
        dfr.gg <- data.frame("grid_index"=gg , "moderator" = moderator.grid[gg] ,
                          "par"= pars0 , "parindex" = 1:NP , dfr.gg    )
        dfr.gg0 <- data.frame("grid_index"=gg , "moderator" = moderator.grid[gg] ,
                          "par"= fit_measures , "parindex" = NP + 1:NF , 
                          "est"= lavaan::fitMeasures(survey.fit , fit.measures= fit_measures ) ,
                          "op"="fit" )
        vars <- setdiff( colnames(dfr.gg) , colnames(dfr.gg0) )
        for (vv in vars){ dfr.gg0[,vv] <- NA }
        dfr.gg <- rbind( dfr.gg , dfr.gg0[ , colnames(dfr.gg) ] )        
        parameters <- rbind( parameters , dfr.gg ) 
        # fits <- rbind( fits , dfr.gg ) 
        if (verbose){
            if ( gg %in% pr ){
                cat("-")
                utils::flush.console()
            }
        }
    }
    if (verbose){
        cat("|\n")
        utils::flush.console()
    }
    parameters <- parameters[ order(parameters$parindex) , ]    
    #--- OUTPUT            
    res <- list( parameters = parameters )
    return(res)    
}
#######################################################################            

lsem.fitsem <- lsem_fitsem
