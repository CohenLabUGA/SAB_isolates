

#### check list of required packages is installed, if not installs and loads them
### function modified from following stackoverflow:
# https://stackoverflow.com/questions/15155814/check-if-r-package-is-installed-then-load-library

checkAndLoadPackages <- function(...,silent=FALSE){
  
  #check names and run 'require' function over if the given package is installed
  requirePkg<- function(pkg){if(length(setdiff(pkg,rownames(installed.packages())))==0)
    require(pkg, quietly = TRUE,character.only = TRUE)
  }
  
  packages <- as.vector(unlist(list(...)))
  if(!is.character(packages))stop("No numeric allowed! Input must contain package names to install and load")
  
  if (length(setdiff(packages,rownames(installed.packages()))) > 0 )
    install.packages(setdiff(packages,rownames(installed.packages())),
                     repos="https://cloud.r-project.org")
  
  res<- unlist(sapply(packages, requirePkg))
  
  if(silent == FALSE && !is.null(res)) {cat("\nBellow Packages Successfully Loaded:\n\n")
    print(res)
  }
}
