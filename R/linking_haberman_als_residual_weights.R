## File Name: linking_haberman_als_residual_weights.R
## File Version: 0.06
## File Last Change: 2017-09-19 21:33:02


linking_haberman_als_residual_weights <- function( logaj , logaAt ,
		logaM , cutoff , wgtM0 , eps )
{
	loga_expected <- TAM::tam_outer( logaj , logaAt , op = "+" )
	loga_resid <- logaM - loga_expected
	wgt_adj <- 1 + 0 * wgtM0
	wgtM <- wgtM0
	if ( cutoff < Inf ){
		wgt_adj <- ( 1 - ( loga_resid / cutoff )^2 )^2
		wgt_adj <- 1 * ( abs( loga_resid ) <= cutoff ) * wgt_adj
		wgtM <- wgtM0 * wgt_adj + eps
	}
	res <- list(loga_resid = loga_resid, wgt_adj = wgt_adj , wgtM = wgtM )
	return(res)		
}
		
