.onAttach <- function(...)
{
  msg <- paste("This is demography", utils::packageVersion("demography"),"\n")
  packageStartupMessage(msg)
}
