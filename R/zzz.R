print.msg.version <- function()
{ library(help=msg)$info[[1]] -> version
  version <- version[pmatch("Version",version)]
  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]
  cat(paste("This is msg",version,"\n"))
}

.onAttach <- function(...) {
  print.msg.version()
}

.onLoad <- function(lib,pkg){
  library.dynam("msg", pkg, lib)
}

.First.lib <- function(lib, pkg) {
  print.msg.version()
}
