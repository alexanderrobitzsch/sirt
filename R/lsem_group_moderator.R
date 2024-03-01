## File Name: lsem_group_moderator.R
## File Version: 0.178


#***** grouping a moderator variable
lsem_group_moderator <- function( data, type, moderator.grid,
            moderator, residualize, h, is_imputed=FALSE, Nimp=0 )
{
    data0 <- data
    moderator.grouped <- NULL
    if (type=='MGM'){
        G1 <- length(moderator.grid)
        moderator.grouped <- data.frame( min=moderator.grid[-G1],
                                            max=moderator.grid[-1] )
        moderator.grouped$mid <- rowMeans( moderator.grouped)
        if (! is_imputed){
            v1 <- data[, moderator ]
            v2 <- moderator.grouped$mid[1]
            for (gg in 2L:G1){
                v2 <- ifelse( v1 > moderator.grouped$max[gg-1],
                            moderator.grouped$mid[gg], v2 )
            }
            data[,moderator] <- v2
        }
        if (is_imputed){
            for (ii in 1L:Nimp){
                data <- data0[[ii]]
                v1 <- data[, moderator ]
                v2 <- moderator.grouped$mid[1]
                for (gg in 2L:G1){
                    v2 <- ifelse( v1 > moderator.grouped$max[gg-1],
                                moderator.grouped$mid[gg], v2 )
                }
                data[,moderator] <- v2
                data0[[ii]] <- data
            }
        }
        # residualize <- FALSE
        h <- 1E-5
        moderator.grid <- moderator.grouped$mid
    }
    res <- list( data=data, moderator.grouped=moderator.grouped,
                residualize=residualize, h=h,
                moderator.grid=moderator.grid )
    return(res)
}


lsem.group.moderator <- lsem_group_moderator
