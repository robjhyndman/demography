.onAttach <- function(...)
{
    library(help=demography)$info[[1]] -> version
    if (!is.null(version)){
        version <- version[pmatch("Version",version)]
        um <- strsplit(version," ")[[1]]
        version <- um[nchar(um)>0][2]
        packageStartupMessage(paste("This is demography",version,"\n"))      
    }
}
