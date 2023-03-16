## File Name: dmlavaan_joint_parameterTable_merge_table.R
## File Version: 0.04

dmlavaan_joint_parameterTable_merge_table <- function(partable, mod, model_index=1)
{
    y1 <- mod$partable[, c('parname','est','se_sw')]
    colnames(y1) <- c('parname',paste0(c('est','se'), model_index ) )
    y1 <- y1[ ! duplicated(y1$parname), ]
    partable <- merge(x=partable, y=y1, by='parname', all.x=TRUE)
    partable <- partable[ order(partable$id), ]
    return(partable)
}
