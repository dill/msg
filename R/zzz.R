.onAttach<-function(...){
  # uses packageStartupMessage which can then be
  # surpressed

  version <- utils::packageVersion("msg")
  built <- utils::packageDescription("msg",fields="Built")

  hello <- paste("This is msg ",version,"\nBuilt: ",built,sep="")
  packageStartupMessage(hello)
}

#.onLoad <- function(lib,pkg){
#  library.dynam("msg", pkg, lib)
#}

