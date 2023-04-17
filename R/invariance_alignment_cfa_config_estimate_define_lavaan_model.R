## File Name: invariance_alignment_cfa_config_estimate_define_lavaan_model.R
## File Version: 0.091


invariance_alignment_cfa_config_estimate_define_lavaan_model <-
        function(items_gg, label_F="F", model="2PM")
{
    while (label_F %in% items_gg){
        label_F <- paste0( label_F, 'F')
    }
    if (substring(model,1,1)=='2'){
        lavmodel <- paste0(label_F, '=~', paste0(items_gg, collapse='+') )
    } else {
        lavmodel <- paste0(label_F, '=~', paste0( paste0('a*',items_gg), collapse='+') )
        lavmodel1 <- paste0(paste0(items_gg, '~~b*', items_gg), collapse='\n' )
        lavmodel <- paste0( lavmodel, '\n', lavmodel1)
    }
    return(lavmodel)
}
