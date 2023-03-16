## File Name: summary.lsem.R
## File Version: 0.408


#-- summary lsem
summary.lsem <- function( object, file=NULL, digits=3, ... )
{
    # open sink for a file
    sirt_osink( file=file  )

    cat("-----------------------------------------------------------------\n")
    cat("Local Structural Equation Model \n\n")

    #-- print packages
    packages <- c("sirt", "lavaan")
    if (object$use_lavaan_survey){
        packages <- c(packages, "lavaan.survey")
    }
    sirt_summary_print_packages(packages=packages)

    #-- print R session
    cat("\n")
    sirt_summary_print_rsession()

    cat(paste0("Function 'sirt::lsem.estimate', type='", object$type,"'"), "\n\n")

    #- print call
    sirt_summary_print_call(CALL=object$CALL)

    #-- print computation time
    sirt_summary_print_computation_time_s1(object=object)

    # space between equality sign
    sp_eq <- paste0( c(" ", "=", " "), collapse="")

    cat( paste0( "Number of observations in datasets", sp_eq,
                    round(object$N, digits) ), "\n")
    cat( paste0( "Used observations in analysis", sp_eq,
                    round(object$nobs, digits) ), "\n")
    cat("Used sampling weights:", ! object$no_sampling_weights, "\n")
    if ( object$type=="LSEM"){
        cat( paste0( "Bandwidth factor", sp_eq, round(object$h,digits) ), "\n")
        cat( paste0( "Bandwidth", sp_eq, round(object$bw,digits) ), "\n")
        cat( paste0( "Number of focal points for moderator", sp_eq,
                            length(object$moderator.grid ) ), "\n")
        cat("\n")
        cat("Used joint estimation:", object$est_joint, "\n")
        cat("Used sufficient statistics:", object$sufficient_statistics, "\n")
        cat("Used local linear smoothing:", object$loc_linear_smooth, "\n")
        cat("Used pseudo weights:", object$use_pseudo_weights, "\n")
        cat("Used lavaan package:", TRUE, "\n")
        cat("Used lavaan.survey package:", object$use_lavaan_survey, "\n\n")
        cat("Mean structure modelled:", object$is_meanstructure, "\n")

        if (object$class_boot){
            v1 <- paste0("\nStatistical inference based on ", object$R,
                                " bootstrap samples.")
            cat(v1,"\n")
        }
    }

    if ( object$type=="MGM"){
        cat( paste0( "Number of groups for moderator=",
                            length(object$moderator.grid ) ), "\n")
    }

    cat("\nlavaan Model\n")
    cat(object$lavmodel)

    if (object$est_joint){
        cat("\n\n")
        cat("Global Fit Statistics for Joint Estimation\n\n")
        obji <- object$fitstats_joint
        sirt_summary_print_objects(obji=obji, digits=digits)
    }

    cat("\n\n")
    cat("Parameter Estimate Summary\n\n")
    obji <- object$parameters_summary
    sirt_summary_print_objects(obji=obji, digits=digits, from=2)

    cat("\n")
    cat("Distribution of Moderator: Density and Effective Sample Size\n\n")
    cat( paste0("M=", round(object$m.moderator, digits), " | SD=",
                round(object$sd.moderator, digits), "\n\n") )
    obji <- object$moderator.density
    sirt_summary_print_objects(obji=obji, digits=digits, from=1)

    cat("\n")
    obji <- object$moderator.stat
    sirt_summary_print_objects(obji=obji, digits=digits, from=2)

    # close file
    sirt_csink(file)
}
