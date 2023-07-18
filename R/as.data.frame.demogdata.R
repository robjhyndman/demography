#' Coerce a demogdata object to a data.frame object
#'
#' @param x Object to be coerced to a data frame.
#' @param ... Other arguments not used
#'
#' @return A data.frame object.
#'
#' @examples
#' # coerce demogdata object to data.frame ----
#' as.data.frame(fr.mort)
#' @export
as.data.frame.demogdata <- function(x, ...) {
  rates_included <- ("rate" %in% names(x))
  pop_included <- ("pop" %in% names(x))
  # Size of matrices
  nyears <- length(x$year)
  nages <- length(x$age)
  if(rates_included)
    groups <- names(x$rate)
  else if(pop_included)
    groups <- names(x$pop)
  else
    groups <- NULL
  outlist <- vector(length=nyears, mode="list")
  # Create data frame for rates
  for(i in seq_along(groups)) {
    outlist[[i]] <- data.frame(
      Year = rep(x$year, rep(nages, nyears)),
      Age = rep(x$age, nyears),
      Group = groups[i]
    )
    if (rates_included) {
      outlist[[i]]$Rates <- c(x$rate[[i]])
    }
    if (pop_included) {
      outlist[[i]]$Exposure <- c(x$pop[[i]])
    }
  }
  out <- do.call("rbind", outlist)
  out$Age <- as.integer(out$Age)
  out$Year <- as.integer(out$Year)
  # Assume Inf rates are due to 0/0
  out$Rates[out$Rates==Inf] <- NA_real_
  # Rename rates column
  if (x$type == "mortality") {
    colnames(out)[4] <- "Mortality"
  } else if (x$type == "fertility") {
    colnames(out)[4] <- "Fertility"
  } else if (x$type == "migration") {
    colnames(out)[4] <- "NetMigration"
  } else {
    stop("Unknown type")
  }
  # Add counts if available
  if (rates_included & pop_included) {
    if ("Mortality" %in% colnames(out) & "Exposure" %in% colnames(out)) {
      out$Deaths <- out$Exposure * out$Mortality
      out$Deaths[is.na(out$Mortality) & out$Exposure > 0] <- 0
    } else if ("Fertility" %in% colnames(out) & "Exposure" %in% colnames(out)) {
      out$Births <- out$Exposure * out$Fertility / 1000
      out$Births[is.na(out$Fertility) & out$Exposure > 0] <- 0
    }
  }
  # Reorganize
  out <- out[order(out$Group, out$Year, out$Age),]
  rownames(out) <- NULL
  return(out)
}

utils::globalVariables(c("Deaths","Births","Year","Age","Exposure","Group","Mortality","Fertility"))

