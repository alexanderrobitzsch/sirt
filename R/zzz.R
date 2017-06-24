#  zzz.R
#
# This function is simply copied from the mice package.


# on attach sirt
.onAttach <- function(libname,pkgname){
  d <- utils::packageDescription("sirt")
  d1 <- d$Version 
#  nk <- paste( rep( " " , 20 - nchar(d1) ) , collapse="")
  packageStartupMessage(
		paste("- " , d$Package," " , d1 ," (",d$Date,")",sep="")  )
	}
	
	
version <- function(pkg="sirt"){
  lib <- dirname( system.file(package = pkg))
  d <- utils::packageDescription(pkg)
  return( paste(d$Package,d$Version,d$Date,lib))
}

# .First.lib <- function(lib, pkg){
#          library.dynam("sirt", package = pkg, lib.loc = lib)
#          return(invisible(0))
#        } 
