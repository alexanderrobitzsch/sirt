## File Name: mirt.wrapper.coef.R
## File Version: 3.154

mirt.wrapper.coef <- function(mirt.obj)
{
    TAM::require_namespace_msg("mirt")
    tmod <- mirt::mod2values(mirt.obj)
    items <- colnames(mirt.obj@Data$data)
    I <- length(items)
    itempars <- NULL
    pars <- list()
    for (ii in 1:I) {
        tmod_ii <- tmod[ tmod$item==items[ii], ]
        pars[[ii]] <- tmod_ii
        itempars <- unique(c(itempars, tmod_ii$name ))
        itempars <- sort(itempars)
        IP <- length(itempars)
    }
    #-- single group case
    if (inherits(mirt.obj,"SingleGroupClass")){
        dfr <- matrix(NA, nrow=I, ncol=IP)
        colnames(dfr) <- itempars
        dfr <- as.data.frame(dfr)
        dfr <- data.frame(item=items, dfr)
        rownames(dfr) <- items
        for (ii in 1:I){
            pars_ii <- pars[[ii]]
            cii <- pars_ii$value
            names(cii) <- pars_ii$name
            dfr[ii, names(cii)] <- cii
        }
        rownames(dfr) <- NULL
        tmod_group <- tmod[ tmod$item=="GROUP", ]
        h2 <- tmod_group$value
        names(h2) <- tmod_group$name
        res <- list(coef=dfr, GroupPars=h2, G=1)
    }
    #-- multiple group case
    if (inherits(mirt.obj,"MultipleGroupClass") ) {
        groups <- sort(unique(tmod$group))
        G <- length(groups)
        dfr0 <- NULL

        dfr <- matrix(NA, nrow=1, ncol=IP)
        colnames(dfr) <- itempars
        dfr <- as.data.frame(dfr)
        dfr0 <- NULL
        for (ii in 1:I){
            for (group in groups){
                pars_ii <- pars[[ii]]
                pars_ii <- pars_ii[ pars_ii$group==group, ]
                cii <- pars_ii$value
                names(cii) <- pars_ii$name
                if (length(cii)>0){
                    dfr1 <- data.frame(group=group, item=items[ii], dfr )
                    dfr1[1, names(cii)] <- cii
                    dfr0 <- rbind(dfr0, dfr1)
                }
            }
        }
        dfr <- dfr0
        dfr <- dfr[ order(dfr$item), ]
        rownames(dfr) <- NULL

        tmod2 <- tmod[ tmod$item=="GROUP", ]
        h1 <- list()
        for (group in groups){
            tmod_group <- tmod2[ tmod2$group==group, ]
            h2 <- tmod_group$value
            names(h2) <- tmod_group$name
            h1[[group]] <- h2
        }
        res <- list(coef=dfr, GroupPars=h1, G=G, groups=groups)
    }
    res$class <- class(mirt.obj)
    return(res)
}
