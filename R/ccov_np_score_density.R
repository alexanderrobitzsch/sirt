## File Name: ccov_np_score_density.R
## File Version: 0.04

ccov_np_score_density <- function(score, thetagrid, smooth=TRUE)
{
    if (smooth){
        thg_dens <- stats::density(x=score, from=min(thetagrid), to=max(thetagrid))
        wgt_thetagrid <- ccov_np_regression(x=thg_dens$x, y=thg_dens$y,
                        xgrid=thetagrid, bwscale=.1)
    } else {
        a1 <- stats::aggregate(1+0*score, list(score), sum, na.rm=TRUE)
        wgt_thetagrid <- a1[,2] / sum(a1[,2])
    }
    return(wgt_thetagrid)
}
