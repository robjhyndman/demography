# Function to construct a mortality demogdata object from HMD
hmd.mx <- function(country, username, password, label=country)
{
    path <- paste("https://www.mortality.org/hmd/", country, "/STATS/", "Mx_1x1.txt", sep = "")
    userpwd <- paste(username, ":", password, sep = "")
    txt <- RCurl::getURL(path, userpwd = userpwd)
    con <- textConnection(txt)
    mx <- try(utils::read.table(con, skip = 2, header = TRUE, na.strings = "."),TRUE)
    close(con)
    if(class(mx)=="try-error")
        stop("Connection error at www.mortality.org. Please check username, password and country label.")

    path <- paste("https://www.mortality.org/hmd/", country, "/STATS/", "Exposures_1x1.txt", sep = "")
    userpwd <- paste(username, ":", password, sep = "")
    txt <- RCurl::getURL(path, userpwd = userpwd)
    con <- textConnection(txt)
    pop <- try(utils::read.table(con, skip = 2, header = TRUE, na.strings = "."),TRUE)
    close(con)
    if(class(pop)=="try-error")
        stop("Exposures file not found at www.mortality.org")

    obj <- list(type="mortality",label=label,lambda=0)

    obj$year <- sort(unique(mx[, 1]))
    #obj$year <- ts(obj$year, start=min(obj$year))
    n <- length(obj$year)
    m <- length(unique(mx[, 2]))
    obj$age <- mx[1:m, 2]
    mnames <- names(mx)[-c(1, 2)]
    n.mort <- length(mnames)
    obj$rate <- obj$pop <- list()
    for (i in 1:n.mort)
    {
        obj$rate[[i]] <- matrix(mx[, i + 2], nrow = m, ncol = n)
        obj$rate[[i]][obj$rate[[i]] < 0] <- NA
        obj$pop[[i]] <- matrix(pop[, i + 2], nrow = m, ncol = n)
        obj$pop[[i]][obj$pop[[i]] < 0] <- NA
        dimnames(obj$rate[[i]]) <- dimnames(obj$pop[[i]]) <- list(obj$age, obj$year)
    }
    names(obj$pop) = names(obj$rate) <- tolower(mnames)

    obj$age <- as.numeric(as.character(obj$age))
    if (is.na(obj$age[m]))
        obj$age[m] <- 2 * obj$age[m - 1] - obj$age[m - 2]
    return(structure(obj, class = "demogdata"))
}


hmd.e0 <- function(country, username, password)
{
    path <- paste("https://www.mortality.org/hmd/", country, "/STATS/", "E0per.txt", sep = "")
    userpwd <- paste(username, ":", password, sep = "")
    txt <- RCurl::getURL(path, userpwd = userpwd)
    con <- textConnection(txt)
    lt <- try(utils::read.table(con, skip = 2, header = TRUE, na.strings = "."), TRUE)
    close(con)
    if (class(lt) == "try-error")
        stop("Life expectancy file not found at www.mortality.org")
	lt <- ts(lt[,-1],start=lt[1,1],frequency=1)
    return(lt)
}


hmd.pop <- function(country, username, password, label=country)
{
    path <- paste("https://www.mortality.org/hmd/", country, "/STATS/", "Population.txt", sep = "")
    userpwd <- paste(username, ":", password, sep = "")
    txt <- RCurl::getURL(path, userpwd = userpwd)
    con <- textConnection(txt)
    pop <- try(utils::read.table(con, skip = 2, header = TRUE, na.strings = "."),TRUE)
    close(con)
    if(class(pop)=="try-error")
        stop("Population file not found at www.mortality.org")

    obj <- list(type="population",label=label,lambda=0)

    obj$year = sort(unique(pop[, 1]))
    #obj$year <- ts(obj$year, start=min(obj$year))
    n <- length(obj$year)
    m <- length(unique(pop[, 2]))
    obj$age <- pop[1:m, 2]
    mnames <- names(pop)[-c(1, 2)]
    n.pop <- length(mnames)
    obj$pop <- list()
    for (i in 1:n.pop)
    {
        obj$pop[[i]] <- matrix(pop[, i + 2], nrow = m, ncol = n)
        obj$pop[[i]][obj$pop[[i]] < 0] <- NA
        dimnames(obj$pop[[i]]) <- list(obj$age, obj$year)
    }
    names(obj$pop) <- tolower(mnames)

    obj$age <- as.numeric(as.character(obj$age))
    if (is.na(obj$age[m]))
        obj$age[m] <- 2 * obj$age[m - 1] - obj$age[m - 2]
    return(structure(obj, class = "demogdata"))
}

