## File Name: xxirt_compute_posterior.R
## File Version: 0.14
## File Last Change: 2017-06-14 19:49:27


###########################################################################
xxirt_compute_posterior <- function( prior_Theta , p.xi.aj , group ,
                 G , weights , dat1 , dat_resp , maxK , group_index,
				 dat1_resp ){
	N <- nrow(dat_resp)
	TP <- ncol(p.xi.aj)	
	I <- ncol(dat1)		
	# posterior distribution
	prior1 <- t( prior_Theta[ , group ] )	 
	p1 <- p.aj.xi <- prior1 * p.xi.aj
	p.aj.xi <- p.aj.xi / rowSums( p.aj.xi )
	# expected counts
	n.ik <- array( 0 , dim=c(I,maxK , TP,G) )		
	N.ik <- array( 0 , dim=c(I,maxK , TP) )		
	pi.k <- matrix( 0 , nrow=TP , ncol=G )		
	for (gg in 1:G){
		ind_gg <- group_index[[gg]]	
		p.aj.xi.gg <- as.matrix( p.aj.xi[ind_gg , ] )
		dat1_resp_gg <- dat1_resp[ ind_gg , , ]
		for (kk in 1:maxK){
			n.ik[,kk,,gg] <- xxirt_compute_posterior_expected_counts(
								dat1_resp_gg = as.matrix(dat1_resp_gg[,,kk] ) ,
								p_aj_xi_gg = p.aj.xi.gg )
		}	
		N.ik <- N.ik + n.ik[,,,gg]
		pi.k[,gg] <- colSums( p.aj.xi.gg * weights[ ind_gg ] )
	}  # end gg
	res <- list( p.aj.xi = p.aj.xi , n.ik = n.ik , N.ik = N.ik, N.k = pi.k ,
					post_unnorm = p1 )
	return(res)
}
###########################################################################			
