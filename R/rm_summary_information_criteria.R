## File Name: rm_summary_information_criteria.R
## File Version: 0.04

rm_summary_information_criteria <- function(object, digits_crit=0, digits_penalty=2)
{        
    ic <- object$ic
    rm_summary_information_criteria_print_one_criterium(ic=ic, crit="AIC", 
            desc="AIC = -2*LL + 2*p ", digits_crit=digits_crit, digits_penalty=digits_penalty)            
    rm_summary_information_criteria_print_one_criterium(ic=ic, crit="AICc", 
            desc="AICc = -2*LL + 2*p + 2*p*(p+1)/(n-p-1)  (bias corrected AIC)", digits_crit=digits_crit, digits_penalty=digits_penalty)            
    rm_summary_information_criteria_print_one_criterium(ic=ic, crit="BIC", 
            desc="BIC = -2*LL + log(n)*p ", digits_crit=digits_crit, digits_penalty=digits_penalty)            
    rm_summary_information_criteria_print_one_criterium(ic=ic, crit="CAIC", 
            desc="CAIC = -2*LL + [log(n)+1]*p  (consistent AIC)", digits_crit=digits_crit, digits_penalty=digits_penalty)            
    cat("\n")
}
