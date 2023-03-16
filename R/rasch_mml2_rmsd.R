## File Name: rasch_mml2_rmsd.R
## File Version: 0.02

rasch_mml2_rmsd <- function(n.jk, r.jk, rprobs, pi.k, dat)
{
    requireNamespace("miceadds")
    dims <- dim(n.jk)
    n.ik <- array(0, dim=c(dims[1],2,dims[2],dims[3]))
    n.ik[,2,,] <- r.jk
    n.ik[,1,,] <- n.jk - r.jk
    n.ik <- aperm(n.ik, perm=c(3,1,2,4))
    probs <- aperm(rprobs, perm=c(3,1,2))
    dim(probs)[[2]] <- ncol(dat)
    rmsd <- CDM::IRT_RMSD_calc_rmsd( n.ik=n.ik, pi.k=pi.k, probs=probs, eps=1E-30 )
    args <- list( n.ik=n.ik, pi.k=pi.k, probs=probs, eps=1E-30 )
    miceadds::Reval( paste0("rmsd$MD <- do.call( CDM:::","IRT_RMSD_calc_md, args=args)"),
                        print.string=FALSE, n.eval.parent=1)
    items <- colnames(dat)
    for (vv in c("RMSD","RMSD_bc","MD","MAD")){
        names(rmsd[[vv]]) <- items
    }
    return(rmsd)
}
