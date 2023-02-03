# demography 2.0

 - Updated hmd functions due to changes at mortality.org. Now using HMDHFDplus for downloads.
 - Added functions to convert lifetable and demogdata objects to data frames.
 - Added pkgdown site

# demography 1.22

 - Made compatible with latest rainbow and ftsa packages

# demography 1.21

 - Using https for HMD
 - roxygenized all documentation
 - made compatible with latest forecast package

# demography 1.20

 - Removed dependency on ftsa now that we no longer need a special median function.

# demography 1.19

 - Lots of clean up to conform to CRAN policy
 - Fixed conflicts with some packages
 - total life expectancy for coherentfdm added
 - Added PI for coherent total life expectancy

# demography 1.18

- Updated lca() with "fertility" data
- Modified handling of warnings in forecast.fdm()
- Fixed problem in simulate() when there are too many missing values in residuals. Now
 all missing residuals are set to 0
- Better handling of weights in forecast.fdm()
- Allowed lca() to handle data that is observed less frequently than annually.

# demography 1.17

- fix for smooth.demogdata caused by changes in mgcv package

# demography 1.16

- Fixed bug in pop.sim when no migration data was used
- Added hmd.pop() to read population data from www.mortality.org.
- Fixed a bug in forecast.fdm() when the time series frequency is greater than 1.
- Added scale argument to read.demogdata()
- Fixed problems in forecasting cohort life expectancy
- Improved documentation for read.demogdata()
- Corrected a bug in smooth.demogdata in the default definition of age.grid

# demography 1.15

- smooth.demogdata will no longer return NAs for fertility data. Instead, the fertility rate for the nearest age with positive rate is used.
- Fixed occasional bug in computing life expectancy prediction intervals from coherent fdm model.
- Changed the way missing values are handled at the ends of the age range when smoothing.
- Allowed missing values when using fdm().

# demography 1.14

- minor changes to lifetable calculation.
- replaced cm.spline and cm.splinefun with wrappers to spline and splinefun, now that these include hyman filtering.

# demography 1.13

- Generalized lca with e0 adjustment to allow starting ages other than 0.
- Modified forecast.fdmpr() to allow better control of fractional differencing parameter.
- Modified tfr to be more robust to series names other than "female".
- Fixed bug in production prediction intervals in flife.expectancy from
 an lca object.
- Added Simon Wood and R Core Team as contributing authors.

# demography 1.12

- Removed partial arg matching throughout.
- Added updating methods for fmforecast and fmforecast2 classes.

# demography 1.11

- added warnings option to forecast.fdm().
- fixed rare error in forecast.fdm()
- fixed incorrect label returned by sex.ratio()

# demography 1.10

- show.labels argument dropped from plot.demogdata() as the facility has been dropped by plot.fds() in the rainbow package.
- Fixed a bug in lifetable() when dealing with age groups of more than 1 year.
- Fixed bug in plot.demogdata() when logarithms of zero rates are calculated.

# demography 1.09-1

- Extended flife.expectancy() for use when there is insufficient historical data to compute cohort life expectancy.
- Fixed a bug in flife.expectancy() when type="cohort".

# demography 1.08

- Fixed bug in plot.fmforecast() when plotting coefficients from lca object.
- Fixed several bugs in flife.expectancy() for forecasting cohort life expectancy.
- Fixed bug in lifetable() when type="cohort" and ages of length 1 to give one additional year.

# demography 1.07

- Modified signs of basis functions and coefficients in fdm() to make interpretation easier. This does not affect final forecasts as the signs cancel.
- Fixed bug in forecast.fdm after fitting with weight=TRUE.
- In forecast.fdmpr(): restricted ARFIMA forecasts for coherent models to use data only from the last K years where K can be specified.

# demography 1.06

- Fixed errors in help file for hmd.e0()
- Fixed bug in forecast.fdm after fitting with weight=TRUE.

# demography 1.05

- Modified lifetable for type="cohort" to prevent partial lifetables being produced unless explicitly requested.
- Added hmd.e0() function.

# demography 1.04

- Fixed bugs in the use of weights in fdm() and smooth.demogdata()
- Lifetable functions rewritten to remove bugs and add additional functionality for cohort lifetables.
- Improved speed of PI calculations in e0
- Added flife.expectancy().

# demography 1.03

- improved documentation for hmd.mx()

# demography 1.02

- changed some examples in the help file for bms() to enable the CRAN checks to run faster.

# demography 1.0

- First version on CRAN
- Added summary() functions for fmforecast, fmforecast2, fdmpr and demogdata objects.
- Added e0 prediction intervals for lca objects
- Added model() functions
- Fixed coherentfdm() to allow use with migration data
- Fixed forecast.fmforecast2 to allow use with migration data
- Fixed simulate.fmforecast2 to allow use with migration data
- Updated pop.sim() to take coherent inputs for mortality and migration
- Added simulation of lca objects
- Fixed lots of bugs
- Changed name of hmd() to hmd.mx() to anticipate other hmd.xx functions in the future.

# demography v0.999 (30 July 2010)

- e0 rewritten to allow calculation from coherentfdm results, and to correct the computation of prediction intervals. These are now done using simulations which are much slower than what was done previously, but they are correct (unlike in previous versions). Set the argument PI=TRUE to compute prediction intervals.
- tfr rewritten to correct the computation of prediction intervals. These are now done using simulations which are much slower than what was done previously, but they are correct (unlike in previous versions). Set the argument PI=TRUE to compute prediction intervals.

# demography v0.998 (21 May 2010)

- Bug fixes in coherentfdm and to make hmd visible.

# demography v0.997 (12 May 2010)

- The package now depends on the ftsa and rainbow packages. All duplicate functions have been omitted.
- A new function hmd allows data to be downloaded directly from the Human Mortality Database.
- A new function coherentfdm and an associated forecast method allows coherent forecasting for groups of functional data.
- Some minor bug fixes.

# demography v0.996 (29 March 2010)

- Fixed bug in lca and added warning to lifetable when there are

# demography v0.995 (4 March 2009)

- Fixed bug in pop.sim.

# demography v0.994 (26 February 2009)

- Corrected lifetable calculations to work when the sex is unknown.

# demography v0.993 (4 August 2008)

- Corrected combine.demogdata to produce a ?pop? object when possible, and modified associated help file accordingly.
- Updated the lifetable function to allow five-??year age groups.

# demography v0.992 (3 July 2008)

- Allowed greater flexibility in fitting stationary coefficients to only some components, and using ar for stationary models.
- Changed the default number of terms in an fdm or ftsm model to 6 rather than 3.
- Bug in smooth.demogdata fixed and made compatible with latest version of mgcv package.
- Added color control in plot.ftsm and plot.fdm

# demography v0.991 (12 May 2008)

- Prediction intervals for fdm objects now allowed using structural time series and random walks with drift.
- Changed default forecasting method to ?arima? in forecast.fdm and forecast.ftsm
- Corrected time component of coefficients from forecast.fdm

# demography v0.99

- Bug fixes in smooth.demogdata
- Bug fix in forecast.fdm to allow method=?ets.na? to work again.

# demography v0.98

- Updated documentation to conform to new CRAN rules.
- Added population forecasting functions as described in Hyndman and Booth (2007).
- Modified the internals of smooth.demogdata to take account of changes in the R base and stats packages.
- A few bug fixes.

# demography v0.97

- Many changes to documentation and functions to satisfy CRAN checks.
- Updated all forecasting functions to work with v1.0 of the forecast package. Check help files as some syntax has changed.
- Various bug fixes.

# demography v0.96

- smooth.demogdata slightly modified to give better results.
- Bug fixes in various functions.
- Smooth.demogdata handles age-??grouping better.
- Smooth.demogdata no longer crashes if all age groups have zero population in a year.

# demography v0.95

- Most data have now been taken out of the package (as it is now publicly available). The only data sets in the demography package are fr.mort (French mortality) and aus.fert (Australian fertility).
- The package now handles fertility, mortality and migration data. The basic data class is ?demogdata? and demogdata$type indicates the type of demographic data.
- Many functions have been revised to handle the new data structures. However, I?ve tried to keep the calling syntax the same. If existing code no longer works, check the help files first.
- fmm has been renamed as fdm as it now fits functional demographic models (and not only functional mortality models).
- read.mortality has been renamed as read.demogdata. Similarly, other functions of the form xxx.mortality are now called xxx.demogdata.
- smooth.demogdata now handles various types of smoothing and smoothing constraints. See the help file for details. It tries to do something appropriate depending on the type of data passed. At the moment, it only handles mortality and fertility data.
- I?ve added tfr to compute total fertility rates, and isfe to compute the Integrated Squared Forecast Error for different model orders.
- The documentation has been revised in many places, and several references added.
- All these changes have almost certainly introduced new errors
