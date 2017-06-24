
######################################################################
# auxiliary function for mirt package
# trace function for all items
Probtrace_sirt <- function(items, Theta)
{
	traces <- lapply(items, mirt::probtrace, Theta=Theta)
	ret <- do.call(cbind, traces)
	return(ret)
}