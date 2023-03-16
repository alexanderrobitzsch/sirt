## File Name: mgsem_proc_model.R
## File Version: 0.286
## File Last Change: 2023-02-06

mgsem_proc_model <- function(model, G=G, random_sd=1e-1, technical, N_group,
        prior_list=NULL, pen_type="lasso", fixed_parms=FALSE,
        partable_start=NULL)
{

    dfr <- NULL
    types <- c("ALPHA", "NU", "LAM", "PHI", "PSI")
    symm_types <- c( "PHI", "PSI")
    N <- sum(N_group)

    names_prior_list <- names(prior_list)
    is_B <- mgsem_proc_model_is_B(model=model)
    technical$is_B <- is_B
    if (is_B){
        types <- c(types, "B")
    }

    I <- mgsem_proc_model_extract_dimension(model=model, entry="est",
                            type="LAM", nrow=TRUE)
    D <- mgsem_proc_model_extract_dimension(model=model, entry="est",
                            type="LAM", nrow=FALSE)

    #** process case of single model
    model <- mgsem_proc_model_single_group(model=model)

    #--- loop over groups
    for (gg in 0:G){

        group <- gg
        hh <- gg+1

        #-- include missing entries
        model[[hh]] <- mgsem_proc_model_include_missing_entries(model_hh=model[[hh]],
                            types=types, entries=c("est","index"), I=I, D=D)
        model_hh <- model[[hh]]
        est <- model_hh$est
        index <- model_hh$index



        for (type in types){

            if (type %in% symm_types){
                symm <- TRUE
            } else {
                symm <- FALSE
            }
            symm0 <- symm
            M1 <- est[[type]]
            if (!is.null(M1)){
                M2 <- index[[type]]
                n1 <- nrow(M1)
                n2 <- ncol(M2)
                for (ii in 1:n1){
                    if (symm){
                        hh <- ii
                    } else {
                        hh <- 1
                    }
                    for (jj in hh:n2){
                        if( M2[ii,jj] !=0 ){   # non-missing entry in 'est'
                            dfr1 <- data.frame( type=type, i1=ii, i2=jj,
                                                group=group)
                            dfr1$name <- paste0(dfr1$type, dfr1$i1, dfr1$i2,
                                                "_G", dfr1$group)
                            dfr1$name2 <- paste0(dfr1$type, dfr1$i1, "-", dfr1$i2,
                                                "_G", dfr1$group)
                            symm <- symm0
                            if (ii==jj){
                                symm <- FALSE
                            }
                            dfr1$symm <- symm
                            dfr1$start <- M1[ii,jj]
                            dfr1$index <- M2[ii,jj]
                            dfr1$est <- dfr1$start
                            dfr1$se <- NA
                            if (gg==0){
                                dfr1$N_group <- N
                            } else {
                                dfr1$N_group <- N_group[gg]
                            }

                            #-- check for entries
                            #-- model specifications
                            entries <- c("lower", "upper","prior", "pen_l2",
                                            "pen_lp", "pen_difflp")
                            dfr1 <- mgsem_proc_model_add_specs_all(model=model_hh,
                                            entries=entries, type=type, ii=ii, jj=jj,
                                            dfr1=dfr1, names_prior_list=names_prior_list,
                                            group=group, N_group=N_group,
                                            pen_type=pen_type)
                            dfr1$unique <- 0
                            dfr1$recycle <- 0
                            #- append to previous parameters
                            dfr <- rbind(dfr, dfr1)
                        }
                    }  # end jj (i2)
                }  # end ii  (i1)
            } # !is.null(M1)
        } #end types

    }    # end gg

    if (any(duplicated(dfr$name))){
        dfr$name <- dfr$name2
    }
    dfr$name2 <- NULL

    #** define parameter indices
    res <- mgsem_proc_model_partable_define_index(partable=dfr)
    dfr <- res$partable
    NP <- res$NP
    ND <- res$ND

    #** define lp entries
    res <- mgsem_proc_model_difflp_information(partable=dfr)
    dfr <- res$partable
    difflp_info <- res$difflp_info
    technical$is_pen_difflp <- difflp_info$is_pen_difflp

    #*** specifications in technical
    technical$is_prior <- sum(dfr$prior!="none") > 0
    technical$is_pen_l2 <- sum(dfr$pen_l2>0) > 0
    technical$is_pen_lp <- sum(dfr$pen_lp>0) > 0

    dfr <- as.data.frame(dfr)

    #** coefficient vector of estimated parameters
    if ( ! is.null(partable_start) ){
        dfr$start <- dfr$est <- partable_start
    }
    coef <- mgsem_partable2coef(partable=dfr)

    #** induce some randomness in starting values
    if (random_sd>0){
        coef <- coef + stats::rnorm(NP, sd=random_sd)
        dfr <- mgsem_coef2partable(coef=coef, partable=dfr)
    }

    #* adapt initial values for bounded estimation
    eps1 <- 1e-2
    dfr$est <- ifelse(dfr$est<dfr$lower, dfr$lower + eps1, dfr$est )
    dfr$est <- ifelse(dfr$est>dfr$upper, dfr$upper - eps1, dfr$est )


    #** include parameter indices
    model <- mgsem_partable2model(partable=dfr, model=model, index=TRUE)

    #*** unique parameters
    loop_parms <- (1:ND)[ dfr$unique==1]

    #- rewrite penalty parameters into model matrices
    entries <- c("pen_l2", "pen_lp", "pen_difflp")
    model <- mgsem_proc_model_update_penalties_matrix(partable=dfr,
                    entries=entries, model=model)

    #--- output
    res <- list(model=model, partable=dfr, NP=NP, ND=ND, coef=coef, I=I, D=D,
                    is_B=is_B, technical=technical, types=types,
                    difflp_info=difflp_info, loop_parms=loop_parms)
    return(res)
}
