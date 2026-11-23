## File Name: hh_linking.R
## File Version: 0.103

hh_linking <- function(mirt_mg, prec_weights="none", rob_weights="none", eps=0.001,
                c=NULL, verbose=TRUE)
{

    # prec_weights, rob_weights or eps can be of length two (first: discrimination,
    #  second: intercept)

    # precision weights prec_weights
    # none or raw (item parameters)

    # robust weights rob_weights
    # none, bisquare1 (significance), bisquare2 (data-driven choice of tuning parameter),
    #  L0, Lp (with values of p)

    # c: tuning factor for bisquare regression

    modlist <- mirt_mg

    prec_weights <- extend_two(x=prec_weights)
    rob_weights <- extend_two(x=rob_weights)
    eps <- extend_two(x=eps)
    c <- extend_two(x=c)

    parinfo <- modlist$parinfo
    G <- modlist$G
    I <- modlist$I
    N <- modlist$N
    items <- modlist$items

    #- create design matrix
    res <- hh_linking_create_design_matrix(modlist=modlist)
    modlist$X <- X <- res$X
    modlist$pwinfo <- pwinfo <- res$parinfo
    modlist$P <- P <- res$P
    modlist$Z <- Z <- res$Z

    #- estimate sigma
    ee <- 1
    res_sigma <- hh_linking_regression_model(modlist=modlist, partype='a',
                        prec_weights=prec_weights[ee], rob_weights=rob_weights[ee],
                        eps=eps[ee], c=c[ee], verbose=verbose)

    modlist$parinfo <- parinfo <- res_sigma$parinfo
    modlist$sigma <- sigma <- res_sigma$delta

    #- estimate mu
    ee <- 2
    res_mu <- hh_linking_regression_model(modlist=modlist, partype='d',
                    prec_weights=prec_weights[ee], rob_weights=rob_weights[ee],
                    eps=eps[ee], c=c[ee], verbose=verbose)

    modlist$parinfo <- parinfo <- res_mu$parinfo
    mu <- res_mu$delta

    delta <- c(mu, sigma)
    names(delta) <- c( paste0('mu', 2L:G), paste0('sigma', 2L:G) )

    #-- output
    res <- list( delta=delta, mu=mu, sigma=sigma, X=X, P=P, Z=Z, G=G, I=I, N=N,
                        parinfo=parinfo, pwinfo=pwinfo, res_mu=res_mu,
                        res_sigma=res_sigma, prec_weights=prec_weights,
                        rob_weights=rob_weights)
    class(res) <- 'hh_linking'
    return(res)
}
