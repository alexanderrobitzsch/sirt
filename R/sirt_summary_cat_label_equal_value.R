## File Name: sirt_summary_cat_label_equal_value.R
## File Version: 0.05
## File Last Change: 2018-12-30



sirt_summary_cat_label_equal_value <- function( label, value, label2="", digits=NULL )
{
    res <- sirt_summary_label_equal_value(label=label, value=value,
                label2=label2, digits=digits)
    cat(res)
}
