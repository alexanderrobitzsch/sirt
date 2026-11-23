## File Name: mirt_mg.R
## File Version: 0.055

mirt_mg <- function(resp, group, weights=NULL, ...)
{
    requireNamespace('mirt')
    groups <- sort(unique(group))
    G <- length(groups)
    if (is.null(weights)){
        weights <- rep(1, nrow(resp))
    }

    #- fit mirt models
    modlist <- list()
    for (gg in 1L:G){
        ind_gg <- which(group==groups[gg])
        resp_gg <- resp[ind_gg,]
        resp_gg <- resp_gg[, colMeans(is.na(resp_gg)) < 1 ]
        modlist[[gg]] <- mirt::mirt( data=resp_gg, model=1, itemtype='2PL',
                                survey.weights=weights[ind_gg], SE=TRUE, ...)
    }

    #- processing
    res <- mirt_mg_proc(modlist=modlist)
    res$modlist <- modlist

    #- output
    class(res) <- 'mirt_mg'
    return(res)
}
