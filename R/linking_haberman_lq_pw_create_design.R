## File Name: linking_haberman_lq_pw_create_design.R
## File Version: 0.141

linking_haberman_lq_pw_create_design <- function(y, ind_studies, ind_items, method)
{
    G <- max(ind_studies)
    I <- max(ind_items)

    des <- data.frame(study=ind_studies, item=ind_items)
    dfr <- NULL
    y_long <- NULL
    X <- NULL
    w <- NULL

    for (gg in 1L:(G-1)){
        for (hh in (gg+1):G){
            items_sel <- intersect( des[ des$study==gg, 'item' ],
                                    des[ des$study==hh, 'item' ] )
            if (length(items_sel)>0){
                dfr1 <- data.frame(study1=gg, study2=hh, item=items_sel)
                I1 <- nrow(dfr1)
                X1 <- matrix(0, I1, G-1)
                if (gg>1){
                    X1[,gg-1] <- 1
                }
                X1[,hh-1] <- -1
                for (ii in items_sel){
                    i1 <- which( ind_studies==gg & ind_items==ii )
                    i2 <- which( ind_studies==hh & ind_items==ii )
                    y1 <- y[i1] - y[i2]
                    y_long <- c(y_long, y1)
                }
                X <- rbind(X, X1)
                dfr <- rbind(dfr, dfr1)
            }
        }  # end hh
    }  # end gg

    #- define weights
    w <- rep(1, nrow(X) )
    if (method=='pw2'){
        a1 <- stats::aggregate(1+0*des$item, list(des$item), sum)
        a2 <- stats::aggregate(1+0*dfr$item, list(dfr$item), sum)
        a1 <- a1[ a1[,1] %in% a2[,1], ]
        w0 <- a1[,2] / a2[,2]
        w <- w0[ match( dfr$item, a1[,1]) ]
    }

    #-- output
    res <- list(y=y_long, X=X, w=w, ind_studies=ind_studies, ind_items=ind_items,
                    design=dfr, G=G, I=I)
    return(res)
}
