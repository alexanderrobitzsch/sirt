## File Name: rm_facets_postproc_person.R
## File Version: 0.06

rm_facets_postproc_person <- function( dat2, dat2.resp, procdata, maxK, RR, theta.k, f.qk.yi )
{
    person <- procdata$person.index
    NP <- nrow(person)
    person$score <- rowSums( dat2 * dat2.resp )
    mkrr <- rep( maxK , RR )
    person$maxscore <- rowSums( dat2.resp * sirt_matrix2( mkrr , nrow=NP) )
    person$EAP <- rowSums( f.qk.yi * sirt_matrix2( theta.k , nrow=NP) )
    person$SE.EAP <- sqrt( rowSums( f.qk.yi * sirt_matrix2( theta.k^2 , nrow=NP) ) - ( person$EAP) ^2 )
    EAP.rel <- rm_eap_reliability( EAP=person$EAP, SE_EAP=person$SE.EAP )
    #--- output
    res <- list( person=person, EAP.rel=EAP.rel )
    return(res)
}
