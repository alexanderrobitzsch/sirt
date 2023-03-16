## File Name: noharm_sirt_number_estimated_parameters.R
## File Version: 0.03
## File Last Change: 2018-12-30

noharm_sirt_number_estimated_parameters <- function(I, Fpatt, Ppatt, Psipatt)
{
    Nestpars <- list(total=0)
    Nestpars$total <- Nestpars$thresh <- I
    Nestpars$F <- sum( Fpatt > 0 )
    Nestpars$total <- Nestpars$total + Nestpars$F
    Nestpars$P <- sum( diag(Ppatt)==1 ) / 2 +  sum( Ppatt==1 ) / 2
    Nestpars$total <- Nestpars$total + Nestpars$P
    Nestpars$Psi <- 1/2 * sum( Psipatt==1 )
    Nestpars$total <- Nestpars$total + Nestpars$Psi
    #--- output
    return(Nestpars)
}
