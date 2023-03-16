## File Name: rasch_jml_person_parameters_summary.R
## File Version: 0.04
## File Last Change: 2018-12-30

rasch_jml_person_parameters_summary <- function(x)
{
    res <- data.frame(N=length(x), mean=mean(x),
                        median=stats::quantile(x, probs=.5),
                        sd=stats::sd(x), min=min(x), max=max(x)
                        )
    return(res)
}
