## File Name: mirt.wrapper.itemplot.R
## File Version: 0.08



mirt.wrapper.itemplot <- function( mirt.obj, ask=TRUE, ...)
{
    TAM::require_namespace_msg("mirt")
    I <- ncol( mirt.obj@Data$data )
    for (ii in 1:I){
        main <- paste0("Trace Lines of Item ", colnames( mirt.obj@Data$data )[ii] )
        print( mirt::itemplot(mirt.obj, item=ii, main=main, ...) )
        graphics::par(ask=ask)
    }
}
