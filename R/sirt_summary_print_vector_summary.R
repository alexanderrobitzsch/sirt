## File Name: sirt_summary_print_vector_summary.R
## File Version: 0.03
## File Last Change: 2018-12-30

sirt_summary_print_vector_summary <- function(obji, digits)
{
    obji <- round( summary(obji), digits)
    print(obji)
}
