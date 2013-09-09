netmigration <- function(mort, fert, mfratio = 1.05)
{
  # Basic checks on inputs
  if (class(mort) != "demogdata" | class(fert) != "demogdata") 
    stop("Inputs not demogdata objects")
  if (mort$type != "mortality") 
    stop("mort not mortality data")
  if (fert$type != "fertility") 
    stop("fert not fertility data")
  
  # Find years with both mortality and fertility data
  yrs <- mort$year[sort(match(fert$year, mort$year))]
  yrs <- mort$year[sort(match(yrs, mort$year - 1))] - 1
  startyearpop <- extract.years(mort, c(yrs, max(yrs) + 1))
  mort <- extract.years(mort, yrs)
  fert <- extract.years(fert, yrs)
  n <- length(yrs)
  p <- length(mort$age)
  
  # Splits births by male and female
  B <- colSums(fert$pop$female * fert$rate$female/1000)
  Bf <- B * 1/(1 + mfratio)
  Bm <- B * mfratio/(1 + mfratio)
  
  # Compute non-survival ratios from mortality rates
  nsr.f <- 1 - lifetable(mort, "female", max.age = max(mort$age))$rx
  nsr.m <- 1 - lifetable(mort, "male", max.age = max(mort$age))$rx
  
  # Estimate deaths
  oldest.f <- colSums(as.matrix(startyearpop$pop$female[(p - 1):p, ]))
  oldest.m <- colSums(as.matrix(startyearpop$pop$male[(p - 1):p, ]))
  if (n > 1)
  {
    Df <- nsr.f * rbind(Bf, startyearpop$pop$female[1:(p - 2), -n - 1], oldest.f[-n - 1])
    Dm <- nsr.m * rbind(Bm, startyearpop$pop$male[1:(p - 2), -n - 1], oldest.m[-n - 1])
  } else
  {
    Df <- nsr.f * c(Bf, startyearpop$pop$female[1:(p - 2), -n - 1], oldest.f[-n - 1])
    Dm <- nsr.m * c(Bm, startyearpop$pop$male[1:(p - 2), -n - 1], oldest.m[-n - 1])
  }
  Dm[is.na(Dm) | Dm < 0] <- 0
  Df[is.na(Df) | Df < 0] <- 0
  
  # Compute net migration
  Mf <- Mm <- matrix(NA, nrow = p, ncol = n)
  for (j in 1:n)
  {
    current <- extract.years(startyearpop, years = yrs[j] + 1)
    prev <- extract.years(startyearpop, years = yrs[j])
    Mf[, j] <- current$pop$female - c(Bf[j], prev$pop$female[1:(p - 2), ], oldest.f[j]) + Df[, j]
    Mm[, j] <- current$pop$male - c(Bm[j], prev$pop$male[1:(p - 2), ], oldest.m[j]) + Dm[, j]
  }
  
  # Store migration figures in a demogdata object with same population figures as startyearpop
  mig <- extract.years(startyearpop, years = yrs)
  mig$rate$female <- Mf
  mig$rate$male <- Mm
  mig$rate$total <- Mm + Mf
  mig$pop$total <- mig$pop$male + mig$pop$female
  mig$lambda <- 1
  dimnames(mig$rate$male) <- dimnames(mig$rate$female) <- dimnames(mig$rate$total) <- dimnames(mig$pop$male)
  mig$type <- "migration"
  
  # Return result. Note: rows are actually (B,0), (0,1), etc.
  return(mig)
}

# data must be a demogdata object containing population values for the last year of observation mort
# and mig are fmforecast2 objects and fert is a fmforecast object.  If they are NULL, it is assumed
# all values are zero.

pop.sim <- function(mort, fert = NULL, mig = NULL, firstyearpop, N = 100, mfratio = 1.05, bootstrap = FALSE)
{
  no.mortality <- FALSE  # Not possible to proceed without mort object
  no.fertility <- is.null(fert)
  no.migration <- is.null(mig)
  
  # Basic checks on inputs
  if (!no.mortality)
  {
    if (class(mort) != "fmforecast2") 
      stop("Inputs not fmforecast2 objects")
    if (mort$female$type != "mortality" | mort$male$type != "mortality") 
      stop("mort not based on mortality data")
  }
  if (!no.fertility)
  {
    if (class(fert)[1] != "fmforecast") 
      stop("Inputs not fmforecast objects")
    if (fert$type != "fertility") 
      stop("fert not based on fertility data")
  }
  if (!no.migration)
  {
    if (class(mig) != "fmforecast2") 
      stop("Inputs not fmforecast2 objects")
    if (mig$male$type != "migration" | mig$female$type != "migration") 
      stop("mig not based on migration data")
  }
  
  firstyr <- mort$male$year
  if (!no.fertility) 
    firstyr <- intersect(firstyr, fert$year)
  if (!no.migration) 
    firstyr <- intersect(firstyr, mig$male$year)
  firstyr <- min(firstyr)
  pop <- extract.years(firstyearpop, firstyr)$pop
  
  # Check ages match.  First make them all integers to prevent integer/numeric clashes
  firstyearpop$age <- as.integer(firstyearpop$age)
  mort$male$age <- as.integer(mort$male$age)
  mig$male$age <- as.integer(mig$male$age)
  if (!no.migration)
  {
    if (!identical(firstyearpop$age, mig$male$age)) 
      stop("Please ensure that migration and population data have the same age dimension")
  }
  if (!no.mortality)
  {
    if (!identical(firstyearpop$age, mort$male$age)) 
      stop("Please ensure that mortality and population data have the same age dimension")
  }
  p <- length(firstyearpop$age)
  
  # Simulate all components
  hm <- hf <- h <- Inf
  if (!no.mortality)
  {
    mort.sim <- simulate(mort, nsim = N, bootstrap = bootstrap)
    hm <- length(mort$male$year)
  }
  if (!no.fertility)
  {
    fert.sim <- simulate(fert, nsim = N, bootstrap = bootstrap)
    hf <- length(fert$year)
  }
  if (!no.migration)
  {
    mig.sim <- simulate(mig, nsim = N, bootstrap = bootstrap)
    h <- length(mig$male$year)
    nm <- length(mig$male$model$year)
  }
  h <- min(hm, hf, h)
  
  # Set up storage space
  if (!no.fertility) 
    fage <- is.element(rownames(pop$female), fert$age)
  pop.f <- pop.m <- array(0, c(p, h, N))
  dimnames(pop.f) <- dimnames(pop.m) <- list(mort$female$age, 1:h, 1:N)
  
  advance <- function(x0, x)
  {
    n <- length(x)
    newx <- c(x0, x[1:(n - 2)], x[n - 1] + x[n])
  }
  
  # Simulate N future sample paths of population numbers
  for (i in 1:N)
  {
    # Start with final observed populations
    popf <- round(c(pop$female))
    popm <- round(c(pop$male))
    
    for (j in 1:h)
    {
      # Compute net migration
      if (no.migration)
      {
        netf <- netm <- 0
        Rf <- popf
        Rm <- popm
      } else
      {
        # Simulate net migration
        netf <- round(mig.sim$female[, j, i] + mig$female$model$res$y[, sample(1:nm, 1)])
        netm <- round(mig.sim$male[, j, i] + mig$male$model$res$y[, sample(1:nm, 1)])
        # Add half migrants to current population
        Rf <- pmax(popf + 0.5 * c(netf[2:(p - 1)], 0.5 * netf[p], 0.5 * netf[p]), 0)
        Rm <- pmax(popm + 0.5 * c(netm[2:(p - 1)], 0.5 * netm[p], 0.5 * netm[p]), 0)
      }
      # Survivorship ratios
      # firstyear pop used only to get structure. Data replaced.
      Rt <- firstyearpop  
      Rt$type <- "mortality"
      Rt$lambda <- 0
      Rt$year <- 1
      Rt$rate$female <- matrix(mort.sim$female[, j, i], ncol = 1)
      Rt$rate$male <- matrix(mort.sim$male[, j, i], ncol = 1)
      Rt$pop$female <- matrix(Rf, ncol = 1)
      Rt$pop$male <- matrix(Rm, ncol = 1)
      Rt$rate$total <- Rt$pop$total <- NULL
      colnames(Rt$pop$male) <- colnames(Rt$pop$female) <- colnames(Rt$rate$male) <- colnames(Rt$rate$female) <- "1"
      rownames(Rt$pop$male) <- rownames(Rt$pop$female) <- rownames(Rt$rate$male) <- rownames(Rt$rate$female) <- Rt$age
      nsr.f <- 1 - lifetable(Rt, "female", max.age = max(Rt$age))$rx
      nsr.m <- 1 - lifetable(Rt, "male", max.age = max(Rt$age))$rx
      
      # Simulate deaths
      cohDf <- pmax(nsr.f[c(2:p, p)] * Rf, 0)
      cohDm <- pmax(nsr.m[c(2:p, p)] * Rm, 0)
      Rf2 <- advance(0, Rf - cohDf)  # Ignore deaths to births for now
      Rm2 <- advance(0, Rm - cohDm)
      Ef <- 0.5 * (Rf + Rf2)
      Em <- 0.5 * (Rm + Rm2)
      Df <- rpois(rep(1, p), Ef * mort.sim$female[, j, i])
      Dm <- rpois(rep(1, p), Em * mort.sim$male[, j, i])
      
      # Compute adjusted population ignoring births
      cohDf[2:(p - 2)] <- 0.5 * (Df[2:(p - 2)] + Df[3:(p - 1)])
      cohDm[2:(p - 2)] <- 0.5 * (Dm[2:(p - 2)] + Dm[3:(p - 1)])
      cohDf[p - 1] <- 0.5 * Df[p - 1] + Df[p]
      cohDm[p - 1] <- 0.5 * Dm[p - 1] + Dm[p]
      cohDf[p] <- cohDm[p] <- 0
      
      Rf2 <- advance(0, Rf - cohDf)  # Fix problem with deaths to births later
      Rm2 <- advance(0, Rm - cohDm)
      
      # Simulate births from Poisson distribution
      if (no.fertility) 
        B <- 0 
      else
      {
        lambda <- 0.5 * (Rf[fage] + Rf2[fage]) * fert.sim[, j, i]/1000
        B <- sum(rpois(rep(1, length(fert$age)), lambda))
      }
      # Randomly split births into two sexes
      Bm <- rbinom(1, B, mfratio/(1 + mfratio))
      Bf <- B - Bm
      
      # Infant mortality
      RfB <- pmax(Bf + 0.5 * netf[1], 0)
      RmB <- pmax(Bm + 0.5 * netm[1], 0)
      cohDfB <- pmax(nsr.f[1] * RfB, 0)
      cohDmB <- pmax(nsr.m[1] * RmB, 0)
      Rf20 <- RfB - cohDfB
      Rm20 <- RmB - cohDmB
      Ef0 <- 0.5 * (Rf[1] + Rf20)
      Em0 <- 0.5 * (Rm[1] + Rm20)
      Df0 <- rpois(1, Ef0 * mort.sim$female[1, j, i])
      Dm0 <- rpois(1, Em0 * mort.sim$male[1, j, i])
      f0f <- cohDfB/(Ef0 * mort.sim$female[1, j, i])
      f0m <- cohDmB/(Em0 * mort.sim$male[1, j, i])
      cohDfB <- f0f * Df0
      cohDmB <- f0m * Dm0
      
      # Now we can fix infant deaths
      cohDf[1] <- (1 - f0f) * Df0 + 0.5 * Df[1]
      cohDm[1] <- (1 - f0m) * Dm0 + 0.5 * Dm[1]
      Rf2 <- advance(RfB - cohDfB, Rf - cohDf)
      Rm2 <- advance(RmB - cohDmB, Rm - cohDm)
      
      # Add remaining migrants
      popf <- round(pmax(Rf2 + 0.5 * netf, 0))
      popm <- round(pmax(Rm2 + 0.5 * netm, 0))
      
      # Store final population numbers
      pop.f[, j, i] <- popf
      pop.m[, j, i] <- popm
    }
  }
  # Set up names of years
  dimnames(pop.m)[[2]] <- dimnames(pop.m)[[2]] <- firstyr + (1:h)-1

  # Return all sample paths
  return(list(male = pop.m, female = pop.f))
}


deaths <- function(x)
{
  if (class(x) != "demogdata" | x$type != "mortality") 
    stop("Not a mortality object")
  npop <- length(x$rate)
  out <- list()
  for (i in 1:npop)
  {
    out[[i]] <- round(x$rate[[i]] * x$pop[[i]])
    out[[i]][is.na(out[[i]])] <- 0
  }
  names(out) <- names(x$rate)
  return(out)
}

births <- function(x)
{
  if (class(x) != "demogdata" | x$type != "fertility") 
    stop("Not a fertility object")
  out <- round(x$rate[[1]] * x$pop[[1]]/1000)
  out[is.na(out)] <- 0
  return(out)
} 
