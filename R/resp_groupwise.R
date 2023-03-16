## File Name: resp_groupwise.R
## File Version: 0.11

resp_groupwise <- function(resp, group, items_group)
{
    dat <- resp
    if (is.null(names(items_group))){
        names(items_group) <- colnames(resp)
    }
    NI <- length(items_group)
    N <- nrow(resp)
    for (ii in seq_len(NI)){
        item_ii <- paste(names(items_group)[ii])
        group_ii <- items_group[[ii]]
        NG <- length(group_ii)
        for (gg in seq_len(NG)){
            ind_gg <- which(group==group_ii[gg])
            v1 <- rep(NA, N)
            v1[ ind_gg ] <- resp[ ind_gg, item_ii]
            dat[ ind_gg, item_ii ] <- NA
            name_gg <- paste0(item_ii, "_Gr", group_ii[gg])
            dat <- data.frame(dat, v1)
            colnames(dat)[ ncol(dat) ] <- name_gg
        }
    }
    dat <- dat[, colSums(!is.na(dat)) > 0 ]
    dat <- dat[, sort(colnames(dat)) ]
    return(dat)
}
