## File Name: sirt_summary_print_vector_summary.R
## File Version: 0.01

sirt_summary_print_vector_summary <- function(obji, digits)
{
    obji <- round( summary(obji), digits)
    print(obji)
}
