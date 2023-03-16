## File Name: ccov_np_score_density.R
## File Version: 0.07
## File Last Change: 2020-02-17

ccov_np_score_density <- function(score, thetagrid, smooth=TRUE)
{
    if (smooth){
        thg_dens <- stats::density(x=score, from=min(thetagrid), to=max(thetagrid))
        wgt_thetagrid <- ccov_np_regression(x=thg_dens$x, y=thg_dens$y,
                                xgrid=thetagrid, bwscale=.1)
        wgt_thetagrid <- sirt_sum_norm(x=wgt_thetagrid)
    } else {
        a1 <- stats::aggregate(x=1+0*score, by=list(score), FUN=sum, na.rm=TRUE)
        wgt_thetagrid <- sirt_sum_norm(x=a1[,2])
    }
    return(wgt_thetagrid)
}
