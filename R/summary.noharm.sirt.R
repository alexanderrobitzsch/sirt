## File Name: summary.noharm.sirt.R
## File Version: 1.188
## File Last Change: 2019-01-07


#--- summary function for noharm.sirt
summary.noharm.sirt <- function( object, file=NULL, ...)
{
    # open sink
    sirt_osink( file=file )

    display_string <- sirt_summary_print_display(symbol="-", len=65)
    cat(display_string)

    #--- package and R session
    sirt_summary_print_package_rsession(pack="sirt")

    #- print call
    sirt_summary_print_call(CALL=object$CALL)

    #-- print computation time
    sirt_summary_print_computation_time_s1(object=object)

    cat( paste0( "Function '", class(object), "'\n\n" ) )

    #-- elapsed time
    sirt_summary_print_elapsed_time(object=object)

    #--- model type
    if (object$modtype==2){
        cat( "Multidimensional Exploratory Factor Analysis\n")
    }
    if (object$modtype==3){
        cat( "Multidimensional Confirmatory Factor Analysis\n")
    }

    if (object$modesttype==1){
        cat( "NOHARM approximation\n")
    }
    # if (object$modesttype==2){
    #     cat( "Estimation based on tetrachoric correlations\n\n")
    # }

    sirt_optimizer_summary_print(res=object$res_opt)

    #*********************
    # descriptives
    cat( paste( "Number of Observations: ", object$Nobs, sep=""), "\n" )
    cat( paste( "Number of Items       : ", object$Nitems, sep=""), "\n" )
    cat( paste( "Number of Dimensions  : ", object$dimensions, sep=""), "\n" )
    # cat( paste( "Number of Iterations  : ", object$iter, sep=""), "\n" )
    if (object$modesttype==1){
        cat( paste( "Tanaka Index          : ", round(object$tanaka,object$display.fit), sep=""), "\n" )
        cat( paste( "RMSR                  : ", round(object$rmsr,object$display.fit), sep=""), "\n\n" )
    } else { cat("\n") }

    cat( paste( "Number of Used Item Pairs      : ", object$sumwgtm, sep=""), "\n" )
    cat( paste( "Number of Estimated Parameters : ", object$Nestpars$total, sep=""), "\n" )
    cat( paste( "       # Thresholds            : ", object$Nestpars$thresh, sep=""), "\n" )
    cat( paste( "       # Loadings              : ", object$Nestpars$F, sep=""), "\n" )
    cat( paste( "       # Variances/Covariances : ", object$Nestpars$P, sep=""), "\n" )
    cat( paste( "       # Residual Correlations : ", object$Nestpars$Psi, sep=""), "\n\n" )

    # chi square statistic
    if (object$modtype %in% 2:4){
        if (object$modesttype==1){
            cat("Chi Square Statistic of Gessaroli & De Champlain (1996)\n\n")
            cat( paste( "Chi2                           : ", round(object$chisquare,3), sep=""), "\n" )
            cat( paste( "Degrees of Freedom (df)        : ", object$df, sep=""), "\n" )
            cat( paste( "p(Chi2,df)                     : ", round(object$p.chisquare,3), sep=""), "\n" )
            cat( paste( "Chi2 / df                      : ", round(object$chisquare / object$df,3), sep=""), "\n" )
            cat( paste( "RMSEA                          : ", round(object$rmsea,3), sep=""), "\n\n" )
        }
        cat( paste( "Green-Yang Reliability Omega Total : ", round(object$omega.rel,3),
                        sep=""), "\n\n" )
    }
    if ( object$modtype %in% 2:4){
        #--- factor correlation
        cat( "Factor Covariance Matrix\n")
        print( round( object$factor.cor,3 ) )
        if ( object$modtype %in% 3){
            cat( "\nFactor Correlation Matrix\n")
            print( round( cov2cor(object$factor.cor), 3 ))
        }
    }
    if ( object$modtype %in% 3:4){
        # item parameter
        l1 <- object$loadings.theta
        cat("\nItem Parameters - Latent Trait Model (THETA) Parametrization\n",
                "Loadings, Constants, Asymptotes and Descriptives\n\n")
        m1 <- l1 %*% as.matrix(object$factor.cor) %*% t( l1 )
        v1 <- 1 + diag(m1)
        l1 <- data.frame( l1, final.constant=object$final.constants,
                    lower=object$lower, upper=object$upper,
                    item.variance=round(v1, 3), N=diag(object$N.itempair),
                    p=diag(object$pm) )
        l1 <- round(l1,3)
        print(l1)
        # residual correlation
        if (object$estpars$estPsi==1){
            cat("\nResidual Correlation Matrix\n")
            l1 <- round( object$residcorr, 3 )
            print(l1)
        }
        # factor analysis parametrization
        l1 <- object$loadings
        cat("\nItem Parameters - Common Factor (DELTA) Parametrization\n",
            "Loadings, Thresholds, Uniquenesses and Asymptotes\n\n")
        l1 <- data.frame( l1, threshold=object$thresholds,
                        lower=object$lower, upper=object$upper,
                        uniqueness=object$uniquenesses )
        l1 <- round(l1,3)
        print( l1 )
    }
    if ( object$modtype %in% 2){
        # item parameters
        # l1 <- object$loadings.theta
        cat("\nItem Parameters - Promax Rotated Parameters (THETA)\n",
                "Loadings, Constants, Asymptotes and Descriptives\n\n")
        l1 <- object$promax.theta
        l1 <- data.frame( l1, final.constant=object$final.constants,
                    lower=object$lower, upper=object$upper,
                    N=diag(object$N.itempair), p=diag(object$pm) )
        l1 <- round(l1,3)
        print(l1)
        res <- list( itempars.promax.theta=l1 )
        #****
        cat("\nItem Parameters - Promax Rotated Parameters (DELTA)\n",
                "Loadings, Constants, Asymptotes and Descriptives\n\n")
        l1 <- object$promax
        l1 <- data.frame( l1, thresh=object$threshold,
                    lower=object$lower, upper=object$upper,
                    N=diag(object$N.itempair), p=diag(object$pm) )
        l1 <- round(l1,3)
        print(l1)
    }

    cat("\n--- Parameter table ---\n\n")
    obji <- object$parm_table
    elim <- c("matid", "starts", "est_par")
    obji <- obji[, ! ( colnames(obji) %in% elim ) ]
    digits <- 3
    sirt_summary_print_objects(obji=obji, digits=digits, from=2, rownames_null=TRUE)

    # close sink
    sirt_csink( file=file )
}
