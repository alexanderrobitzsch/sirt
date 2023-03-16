## File Name: mml_raschtype_counts.R
## File Version: 0.08



# calculation of counts
mml_raschtype_counts <- function (dat2, dat2resp, dat1, fqkyi, pik, fyiqk)
{
    res <- MML2_RASCHTYPE_COUNTS( DAT2=dat2, DAT2RESP=dat2resp, DAT1=dat1,
                FQKYI=fqkyi, PIK=pik, FYIQK=fyiqk)
    return(res)
}
