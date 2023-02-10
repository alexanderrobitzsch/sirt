## File Name: lsem_bootstrap_draw_bootstrap_sample.R
## File Version: 0.067

lsem_bootstrap_draw_bootstrap_sample <- function(data, sampling_weights,
    lsem_args, cluster=NULL, repl_design=NULL, rr=NULL)
{
    lsem_args1 <- lsem_args
    ind <- NULL
    used_repl_design <- ! is.null(repl_design)
    repl_vector <- NULL
    N1 <- N <- nrow(data)

    if (! used_repl_design){
        # no replication design
        if (is.null(cluster)){
            ind <- sample(1:N, N, replace=TRUE)
        } else {
            cluster1 <- data[,cluster]
            t1 <- unique(cluster1)
            N <- length(t1)
            ind0 <- sort(sample(t1, size=N, replace=TRUE))
            ind <- NULL
            for (nn in 1L:N){
                v1 <- which(cluster1==ind0[nn])
                ind <- c(ind, v1)
            }
        }

        lsem_args1$data <- data[ind,]
        lsem_args1$sampling_weights <- sampling_weights[ind]

        # define replication design
        t1 <- table(ind)
        t1 <- data.frame( index=as.numeric(names(t1)),
                                freq=as.numeric(t1) )
        rownames(t1) <- NULL

        t2 <- data.frame(index=1:N1, weight=sampling_weights)
        t1 <- merge(x=t1, y=t2, by="index", all=TRUE)
        t1$freq <- ifelse( is.na(t1$freq), 0, t1$freq )
        t1$repl_vector <- t1$freq * t1$weight
        repl_vector <- t1$repl_vector

    } else { # replication design
        lsem_args1$sampling_weights <- repl_design[,rr]
        repl_vector <- lsem_args1$sampling_weights
    }

    #--- arrange output
    res <- list(lsem_args1=lsem_args1, repl_vector=repl_vector,
                    used_repl_design=used_repl_design)

    return(res)
}
