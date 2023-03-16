## File Name: IRT.posterior.mirt.R
## File Version: 0.21
## File Last Change: 2019-10-29



#--- likelihood singleGroup
IRT.posterior.mirt.singleGroup <- function( object, ... )
{
    ll <- IRT_likelihood_mirt_singleGroup( object=object, type="f.qk.yi", ...)
    return(ll)
}
IRT.posterior.ConfirmatoryClass <- IRT.posterior.mirt.singleGroup
IRT.posterior.ExploratoryClass <- IRT.posterior.mirt.singleGroup
IRT.posterior.SingleGroupClass <- IRT.posterior.mirt.singleGroup

IRT.posterior.mirt.multipleGroup <- function( object, ... )
{
    ll <- IRT_likelihood_mirt_multipleGroup( object=object, type="f.qk.yi", ... )
    return(ll)
}
IRT.posterior.MultipleGroupClass <- IRT.posterior.mirt.multipleGroup
