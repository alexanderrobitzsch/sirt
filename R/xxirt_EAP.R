
############################################
# compute EAP and its standard deviation
xxirt_EAP <- function(p.aj.xi , Theta )
{
	D <- ncol(Theta)
	e1 <- p.aj.xi %*% Theta
	colnames(e1) <- paste0("EAP.Dim",1:D)
	e2 <- p.aj.xi %*% Theta^2
	colnames(e2) <- paste0("SD.EAP.Dim",1:D)	
	e2 <- e2 - e1^2 
	res <- cbind( e1 , e2)
	res <- res[ , rep( c(0,D) , D ) + 1:D ]
	return(res)
}
################################################		
