## File Name: rasch_mml2_prior_information_generate_string.R
## File Version: 0.02


rasch_mml2_prior_information_generate_string <- function(prior, distribution)
{
    string <- "None"
    if (!is.null(prior)){
        string <- paste0(distribution,"(", prior[1], ",", prior[2],")")
    }
    return(string)
}
