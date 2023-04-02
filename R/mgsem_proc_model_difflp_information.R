## File Name: mgsem_proc_model_difflp_information.R
## File Version: 0.187


mgsem_proc_model_difflp_information <- function(partable)
{
    ND <- nrow(partable)
    NP <- max(partable$index)
    is_lpdiff_entry <- rep(FALSE,ND)
    lpdiff_diff_indices <- list()
    pen_difflp_factor <- list()

    pen_entry <- 'pen_difflp'
    loop_parms <- which(partable$unique==1 & partable$group>0)

    for (dd in loop_parms){
        index_dd <- partable[dd,'index']
        ind <- which( partable$type==partable$type[dd] & partable$i1==partable$i1[dd] &
                partable$i2==partable$i2[dd] & (partable$group > 0) & (partable$rid!=dd) &
                partable$unique==1)
        if (partable$group[dd]==0){
            ind <- NULL
        }
        if ( (length(ind)>0) & (partable$pen_difflp[dd]>0)  ){
            is_lpdiff_entry[dd] <- TRUE
            lpdiff_diff_indices[[dd]] <- partable[ind, 'index']
            fac <- sqrt(partable[dd,pen_entry]*partable[ind,pen_entry])
            pen_difflp_factor[[dd]] <- fac
        } else {
            partable$pen_difflp[dd] <- 0
        }
    }

    #--- design matrices
    cp <- partable[ partable$unique==1, 'name']
    lpdiff_facmat <- matrix(0, nrow=NP, ncol=NP)
    rownames(lpdiff_facmat) <- colnames(lpdiff_facmat) <- cp
    lpdiff_n <- lpdiff_facmat
    for (dd in loop_parms){
        if (is_lpdiff_entry[dd]){
            i1 <- partable[dd,'index']
            i2 <- lpdiff_diff_indices[[dd]]
            lpdiff_facmat[i1,i2] <- pen_difflp_factor[[dd]]
            n1 <- partable[partable$index==i1,'N_group']
            n2 <- partable[partable$index%in%i2,'N_group']
            lpdiff_n[i1,i2] <- 0.5*sqrt(n1*n2)
        }
    }

    difflp_indices <- which( rowSums(lpdiff_facmat)>0 )
    lpdiff_facmat <- lpdiff_facmat[difflp_indices, difflp_indices]
    lpdiff_n <- lpdiff_n[difflp_indices, difflp_indices]
    is_pen_difflp <- nrow(lpdiff_facmat)>0
    lpdiff_facmat_logical <- (lpdiff_facmat>0)

    #--- output
    difflp_info <- list(is_lpdiff_entry=is_lpdiff_entry,
                        lpdiff_diff_indices=lpdiff_diff_indices,
                        pen_difflp_factor=pen_difflp_factor,
                        lpdiff_facmat=lpdiff_facmat,
                        difflp_indices=difflp_indices,
                        is_pen_difflp=is_pen_difflp,
                        lpdiff_facmat_logical=lpdiff_facmat_logical,
                        lpdiff_n=lpdiff_n )
    res <- list(partable=partable, difflp_info=difflp_info)
    return(res)
}
