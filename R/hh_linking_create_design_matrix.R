## File Name: hh_linking_create_design_matrix.R
## File Version: 0.093


hh_linking_create_design_matrix <- function(modlist)
{
    parinfo <- modlist$parinfo
    G <- modlist$G
    I <- modlist$I
    items <- modlist$items
    parinfo1 <- parinfo[ parinfo$partype=='a', ]
    NP <- nrow(parinfo1)

    X <- matrix(0, nrow=NP, ncol=G-1+I)
    colnames(X) <- c( paste0('G',2L:G), items )
    rownames(X) <- parinfo1$par

    for (gg in 2L:G){
        X[ parinfo1$group==gg, gg-1] <- 1
    }
    for (ii in 1L:I){
        X[ parinfo1$itemid==ii, G-1+ii] <- 1
    }

    #- design matrix for pairwise modeling

    Z <- NULL
    P <- NULL
    dfr <- NULL

    gg <- 1
    hh <- 2

    for (gg in 1L:(G-1)){
        for (hh in (gg+1):G){
            parinfo1_gg <- parinfo1[ parinfo1$group==gg, ]
            parinfo1_hh <- parinfo1[ parinfo1$group==hh, ]
            itemid_gg_hh <- intersect(parinfo1_gg$itemid, parinfo1_hh$itemid)
            I1 <- length(itemid_gg_hh)
            if (I1>0){
                parinfo1_gg <- parinfo1_gg[ parinfo1_gg$itemid %in% itemid_gg_hh,]
                parinfo1_hh <- parinfo1_hh[ parinfo1_hh$itemid %in% itemid_gg_hh,]

                dfr1a <- data.frame(group1=gg, group2=hh, item=parinfo1_gg$item,
                                        itemid=parinfo1_gg$itemid,
                                        parid1=parinfo1_gg$parid,
                                        parid2=parinfo1_hh$parid)
                P1 <- matrix(0, nrow=I1, ncol=NP)
                colnames(P1) <- parinfo1$par
                P1[ cbind(1L:I1, dfr1a$parid1) ] <- 1
                P1[ cbind(1L:I1, dfr1a$parid2) ] <- -1
                Z1 <- matrix(0, nrow=I1, ncol=G-1)
                colnames(Z1) <- paste0('G',2L:G)
                if (gg>1){
                    Z1[, gg-1] <- 1
                }
                Z1[,hh-1] <- -1
                Z <- rbind(Z, Z1)
                P <- rbind(P, P1)
                dfr <- rbind(dfr, dfr1a)
            }
        }
    }

    #-- output
    res <- list(X=X, Z=Z, P=P, pwinfo=dfr)
    return(res)
}
