## File Name: print_digits.R
## File Version: 0.03
## File Last Change: 2023-03-16

print_digits <- function(x, digits=NULL)
{
    NC <- ncol(x)
    if (length(digits)==1){
        digits <- rep(digits, NC)
    }
    for (cc in 1:NC){
        y <- x[,cc]
        if (is.numeric(y)){
            x[,cc] <- round(y, digits=digits[cc])
        }
    }
    print(x)
    invisible(x)
}
