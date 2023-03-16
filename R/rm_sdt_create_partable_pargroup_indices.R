## File Name: rm_sdt_create_partable_pargroup_indices.R
## File Version: 0.08
## File Last Change: 2018-12-30


rm_sdt_create_partable_pargroup_indices <- function(partable, item.index,
    diffindex)
{
    pargroup_diffindex <- list()
    pargroup_index <- list()
    pargroup_type <- list()
    max_pargroup <- max(partable$pargroup)
    for (pp in seq_len(max_pargroup) ){
        pt_pp <- partable[ partable$pargroup==pp, ]
        if (pt_pp$type[1]=="tau"){
            I <- max(partable$row)
            g1 <- seq_len(I)
            item0 <- which( pt_pp$parindex < 0 )
            g1[ g1 %in% item0 ] <- 0
            pargroup_diffindex[[pp]] <- g1
        }
        if (pt_pp$type[1]=="a"){
            g1 <- seq_len(I)
            partable_a <- partable[ partable$type=="a", ]
            if ( sd( partable_a$parindex)==0 ){
                g1 <- 1+0*g1
            }
            item0 <- which( pt_pp$parindex < 0 )
            g1[ g1 %in% item0 ] <- 0
            pargroup_diffindex[[pp]] <- g1
        }
        if (pt_pp$type[1]=="c"){
            g1 <- diffindex$c.rater
            item0 <- which( pt_pp$parindex < 0 )
            g1[ item.index %in% item0 ] <- 0
            pargroup_diffindex[[pp]] <- g1
        }
        if (pt_pp$type[1]=="d"){
            g1 <- diffindex$d.rater
            item0 <- which( pt_pp$parindex < 0 )
            g1[ item.index %in% item0 ] <- 0
            pargroup_diffindex[[pp]] <- g1
        }
        p1 <- partable[ partable$est & partable$pargroup==pp, "parindex" ]
        pargroup_index[[pp]] <- p1
        pargroup_type[[pp]] <- pt_pp$type[1]
    }
    max_pargroup <- max( partable$pargroup )
    np <- max(partable$parindex)

    #--- output
    res <- list( pargroup_diffindex=pargroup_diffindex, np=np,
                pargroup_index=pargroup_index, pargroup_type=pargroup_type,
                max_pargroup=max_pargroup)
    return(res)
}
