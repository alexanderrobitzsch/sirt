## File Name: move_variables_df.R
## File Version: 0.02

##-- move variables in a data frame
move_variables_df <- function(x, after_var, move_vars)
{
    cnx <- colnames(x)
    i1 <- which(cnx==after_var)
    vars1 <- cnx[ seq(1,i1) ]
    vars2 <- move_vars
    vars3 <- setdiff( cnx, c(vars1, vars2) )
    x <- x[, c(vars1, vars2, vars3) ]
    return(x)
}
