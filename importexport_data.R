# This code is licensed under the MIT License. See the LICENSE file for details.


# This file is used for loading and preparing the data for the analysis.



# improve error handling
options(error = function() {
  sink(stderr())
  on.exit(sink(NULL))
  traceback(3, max.lines = 1L)
  if (!interactive()) {
    q(status = 1)
  }
})

# load libraries
library("eurostat")
library("knitr")
library("dplyr")
library("fredr")
library("forecast")
library("xts")
library("tidyverse")
library("marima")
library("MTS")
library("TSA")
library("doSNOW")
library("beepr")

source(".\\importexport_functions.R")

# -----  determine which data should be downloaded
id <- "ei_eteu27_2020_m"
dat <- get_eurostat(id, time_format="num", stringsAsFactors=TRUE)
dat <- data.frame(dat)

ctry <- 'DE'
tradetype <- 'EXP'
mypartner <- 'EU27_2020'
mypartner <- 'EXT_EU27_2020'

infl_start_year <- 1998


# ----- load the data

# country, partner, import specifically
dat.filt <- filter(dat, indic=='ET-T', geo == ctry, unit=='MIO-EUR-NSA', stk_flow==tradetype, partner==mypartner)

y <- ts(dat.filt$values, frequency=12, start=dat.filt$TIME_PERIOD[1], end=dat.filt$TIME_PERIOD[nrow(dat.filt)])
plot(y)

# country, partner specifically
dat.filtctry <- filter(dat, indic=='ET-T', geo == ctry, unit=='MIO-EUR-NSA', partner==mypartner)
y.ctry <- dat.filtctry %>%
  pivot_wider(
    names_from = stk_flow,  
    values_from = values        
  )
y.ctry[, 1:5] <- NULL
y.ctry <- as.data.frame(y.ctry)
head(y.ctry)


# partner and import specifically
mycountries <- c("DE", "FR", "NL", "BE", "DK")
dat.filtexp <- filter(dat, indic=='ET-T', unit=='MIO-EUR-NSA', stk_flow == "EXP", partner==mypartner, geo %in% mycountries)
y.exp <- dat.filtexp %>%
  pivot_wider(
    names_from = geo,  
    values_from = values              
  )


y.exp[, 1:5] <- NULL
y.exp <- as.data.frame(y.exp)
head(y.exp)

# ----- reorganize the data to have one column per country

dat.filtexp.all <- filter(dat, indic=='ET-T', unit=='MIO-EUR-NSA', stk_flow == "EXP", partner=="EXT_EU27_2020")
y.exp.all <- dat.filtexp.all %>%
  pivot_wider(
    names_from = geo,  
    values_from = values              
  )

dat.filtimp.all <- filter(dat, indic=='ET-T', unit=='MIO-EUR-NSA', stk_flow == "IMP", partner=="EXT_EU27_2020")
y.imp.all <- dat.filtimp.all %>%
  pivot_wider(
    names_from = geo,  
    values_from = values              
  )

dat.filtexp.int.all <- filter(dat, indic=='ET-T', unit=='MIO-EUR-NSA', stk_flow == "EXP", partner=="EU27_2020")
y.exp.int.all <- dat.filtexp.int.all %>%
  pivot_wider(
    names_from = geo,  
    values_from = values              
  )

ctrynames <- unique(dat.filtexp.all$geo)
ctrynamesHR <- ctrynames[-11]
ctrynames <- ctrynames[c(-14,-11)] # remove HR and EU2020


# ----- construct training set


# start of training set is first available date
start_train <- min(y.exp.all$TIME_PERIOD)


# ----- principal components for trend prediction

pc.exp.all <- prcomp(y.exp.all[,colnames(y.exp.all) %in% ctrynames], scale. = TRUE)
pc.imp.all <- prcomp(y.imp.all[,colnames(y.imp.all) %in% ctrynames], scale. = TRUE)
pc.exp.int.all <- prcomp(y.exp.int.all[,colnames(y.exp.int.all) %in% ctrynames], scale. = TRUE)


# ----- load inflation data

sysdate <- Sys.Date()
sysyr <- format(sysdate, "%Y")
inflation_filters = list(unit = "I15", freq = "M", coicop = "CP00", untilTimePeriod=sysyr, sinceTimePeriod=infl_start_year)
inflation_data <- get_eurostat(id = "prc_hicp_midx", time_format = "num", stringsAsFactors = TRUE, filters=inflation_filters)
#inflation_filtered <- filter(inflation_data, geo %in% mycountries)
inflation_filtered <- inflation_data

infl.dat <- inflation_filtered %>%
  pivot_wider(
    names_from = geo,  
    values_from = values              
  )

# ----- predict final inflation point
ninfl <- nrow(infl.dat)
infl.oct <- infl.dat[ninfl,]
infl.oct$time <- infl.oct$time+1/12
for (k.ctry in ctrynamesHR) {
  infl.ctry <- infl.dat[k.ctry]
  inflation_ts <- ts(infl.ctry, frequency = 12, start = c(infl_start_year,1))
  infl.oct[k.ctry] <- tail(inflation_ts,1) * geom.mean(quotients(inflation_ts))
}
infl.dat[ninfl+1,] <- infl.oct




# ----- compute inflation-adjusted and pca
y.exp.all.infl <- y.exp.all
y.imp.all.infl <- y.imp.all
y.exp.int.all.infl <- y.exp.int.all
data_start <- min(y.exp.all$TIME_PERIOD)
for (k.ctry in ctrynames) {
  infl.ctry <-  ts(rbind(infl.dat[k.ctry], infl.dat[nrow(infl.dat), k.ctry]),
                   frequency = 12, start = c(infl_start_year,1))
  # create a time series to divide across the correct timestamps and save into tibble column
  y.exp.all.infl[k.ctry] <- as.numeric( ts(y.exp.all[k.ctry], frequency=12, start=data_start) / infl.ctry )
  y.imp.all.infl[k.ctry] <- as.numeric( ts(y.imp.all[k.ctry], frequency=12, start=data_start) / infl.ctry )
  y.exp.int.all.infl[k.ctry] <- as.numeric( ts(y.exp.int.all[k.ctry], frequency=12, start=data_start) / infl.ctry )
}

pc.exp.all.infl <- prcomp(y.exp.all.infl[,colnames(y.exp.all) %in% ctrynames], scale. = TRUE)
pc.imp.all.infl <- prcomp(y.imp.all.infl[,colnames(y.imp.all) %in% ctrynames], scale. = TRUE)
pc.exp.int.all.infl <- prcomp(y.exp.int.all.infl[,colnames(y.exp.int.all) %in% ctrynames], scale. = TRUE)



# ----- create dataframes of used pca values

pca.df <- data.frame(time=y.exp.all$TIME_PERIOD)

pca.df$pc.e.1. <- as.numeric(computeTrend(ts(pc.exp.all.infl$x[,1], frequency=12, start=start_train), 1, 12))
pca.df$pc.ein.1. <- as.numeric(computeTrend(ts(pc.exp.int.all.infl$x[,1], frequency=12, start=start_train), 1, 12))
pca.df$pc.i.1. <- as.numeric(computeTrend(ts(pc.imp.all.infl$x[,1], frequency=12, start=start_train), 1, 12))

pca.df$pc.e.2. <- as.numeric(computeTrend(ts(pc.exp.all.infl$x[,2], frequency=12, start=start_train), 1, 12))
pca.df$pc.ein.2. <- as.numeric(computeTrend(ts(pc.exp.int.all.infl$x[,2], frequency=12, start=start_train), 1, 12))
pca.df$pc.i.2. <- as.numeric(computeTrend(ts(pc.imp.all.infl$x[,2], frequency=12, start=start_train), 1, 12))


# ----- generate pca on stationary values

# generate noise pca: second version: pca on noise
pca.noise.df <- data.frame(time=y.exp.all$TIME_PERIOD)

pca.temp <- compute_noise_pca_df(y.exp.all.infl, ctrynames, 12)
pca.noise.df$pc.e.1. <- c(rep(NA, 11), pca.temp$x[,1])
pca.noise.df$pc.e.2. <- c(rep(NA, 11), pca.temp$x[,2])
#pca.noise.df$pc.e.3. <- c(rep(NA, 11), pca.temp$x[,3])

pca.temp <- compute_noise_pca_df(y.exp.int.all.infl, ctrynames, 12)
pca.noise.df$pc.ein.1. <- c(rep(NA, 11), pca.temp$x[,1])
pca.noise.df$pc.ein.2. <- c(rep(NA, 11), pca.temp$x[,2])
#pca.noise.df$pc.ein.3. <- c(rep(NA, 11), pca.temp$x[,3])

pca.temp <- compute_noise_pca_df(y.imp.all.infl, ctrynames, 12)
pca.noise.df$pc.i.1. <- c(rep(NA, 11), pca.temp$x[,1])
pca.noise.df$pc.i.2. <- c(rep(NA, 11), pca.temp$x[,2])
#pca.noise.df$pc.i.3. <- c(rep(NA, 11), pca.temp$x[,3])


