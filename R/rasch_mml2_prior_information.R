## File Name: rasch_mml2_prior_information.R
## File Version: 0.04
## File Last Change: 2019-10-27

rasch_mml2_prior_information <- function(prior.a, prior.b, prior.c,
    prior.d)
{
    a <- rasch_mml2_prior_information_generate_string(prior=prior.a,
                    distribution="N")
    b <- rasch_mml2_prior_information_generate_string(prior=prior.b,
                    distribution="N")
    c <- rasch_mml2_prior_information_generate_string(prior=prior.c,
                    distribution="Beta")
    d <- rasch_mml2_prior_information_generate_string(prior=prior.d,
                    distribution="Beta")
    priors <- list( a=a, b=b, c=c, d=d)
    return(priors)

}
