.onAttach <- function(...)
{
  msg <- paste("This is demography", packageVersion("demography"),"\n")
  packageStartupMessage(msg)
}
