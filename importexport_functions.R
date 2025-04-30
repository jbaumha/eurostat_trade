# This code is licensed under the MIT License. See the LICENSE file for details.



library("forecast")
library("marima")
library("MTS")
library("TSA")


compute_rmse <- function(y.dat, start_train, start_test, end_test, p = 1, q = 0, d = 0, h = 1, 
                         useinflation = FALSE, trendfc.type = 1, noisefc.type = 1, 
                         inflation.dat = NULL, plotPredictions = FALSE,
                         pca.trend.dat = NULL, pca.noise.dat = NULL) {
  # maybe improve efficiency by only recomputing trend / seasonality for new entry in training data
  
  sqsum <- 0
  nErr <- 0
  
  # dates we evaluate at
  y.alltests <- window(y.dat, start=start_test, end=end_test)
  testdates <- time(y.alltests)
  ntests <- length(testdates)
  
  rses <- y.dat
  rses[] <- 0
  
  if (plotPredictions) {
    y.forecasts <- y.dat
  }
  
  for (testsetDate in testdates) {
    y.train <- window(y.dat, start=start_train, end=testsetDate-h/12)
    y.test <- window(y.dat, start=testsetDate, end=testsetDate)
    
    
    if (useinflation) {
      infl.train <-  window(inflation.dat, start=start_train, end=testsetDate-h/12)
      y.train <- y.train/infl.train
      infl.test <- as.numeric(window(inflation.dat, start=testsetDate-h/12, end=testsetDate-h/12))
      y.test <- y.test/infl.test
    }
    
    if (trendfc.type==-1) { # use sarima
      arima.mod <- arima(y.train, order=c(p,d,q), seasonal=list(order=c(1,0,0), period=12))
      y.fc <- predict(arima.mod,1)$pred
    } else {
    
      # detrend training data
      filter.size <- 12
      y.train.trend <- computeTrend(y.train, 1, filter.size)
      y.train.detrended <- y.train - y.train.trend
      y.train.seasyr <- computeSeas(y.train.detrended)
      y.train.rand <- removeSeas(y.train.detrended, y.train.seasyr)
      
      # fit model and forecast
      if (noisefc.type == 0) {
        fc <- 0
      } else if (noisefc.type == 1) {
        fc <- Inf
        result = tryCatch({
          fc <- forecast(Arima(y.train.rand, order=c(p,0,q)), h = h)$mean[1]
        }, error = function(e) {
          
        })
      } else if (noisefc.type == 2) {
        tryCatch({ # in case system is "exactly singular", forecast 0.
        # use PCA predictors
        
        # build pc.df
        pc.noise.train <- pca.noise.dat[pca.noise.dat$time <= (testsetDate-h/12+1e-8),]
        pc.df <- data.frame(time = pc.noise.train$time)
        pc.noise.train <- pc.noise.train[,-1]
        
        # append all lag 2 from pca
        for (k.colname in colnames(pc.noise.train)) {
          pc.df[paste(k.colname,1,sep="")] <- dplyr::lag(pc.noise.train[k.colname],h)
          pc.df[paste(k.colname,2,sep="")] <- dplyr::lag(pc.noise.train[k.colname],h+1)
        }
        pc.df <- pc.df[,-1]
        
        ax.mod <- arimax(y.train.rand, order=c(p,d,q), xreg=pc.df)
        
        newxreg <- data.frame(time=1)
        for (k.colname in colnames(pc.noise.train)) {
          newxreg[paste(k.colname,1,sep="")] <- tail(pc.noise.train[k.colname],1)
          newxreg[paste(k.colname,2,sep="")] <- tail(pc.noise.train[k.colname],2)[1,]
        }
        newxreg <- newxreg[,-1]
        
        fc <- as.numeric(predict(ax.mod, n.ahead=h, newxreg)$pred)[h]
        }, error = function(e) {
          fc <- Inf
        })
      } else {
        fc <- 0
      }
      
      # if fit or forecast unsuccessful, forecast noise with 0
      if (fc == Inf) {
        fc <- 0
        nErr <- nErr + 1
      }
      
      # compute seasonal component and trend for forecast
      y.fc.month <- round((time(y.test)[1]%%1)*12)
      y.fc.seas <- y.train.seasyr[y.fc.month+1]
      
      if (trendfc.type == 1) {
        y.fc.trend <- tail(y.train.trend, n=1)
      } else if (trendfc.type == 2) {
        y.fc.trend.dat <- window(y.train, start=testsetDate-(h-1+filter.size)/12, end=testsetDate-h/12)
        y.fc.trend.df <- data.frame(xval = time(y.fc.trend.dat), yval = y.fc.trend.dat)
        colnames(y.fc.trend.df) <- c("xval", "yval")
        y.fc.trend.mod <- lm(yval ~poly(xval, 2), data=y.fc.trend.df)
        y.fc.trend <- predict(y.fc.trend.mod, data.frame(xval=time(y.test)[1]))
      } else if (trendfc.type==3) {
        # use last 2 values from trend (after shifting by h)
        y.tr.num <- as.numeric(y.train.trend)
        y.fc.trend.df <- data.frame(yval = y.tr.num, 
                                    yval.lagH = lag(y.tr.num, h), 
                                    yval.lagH1 = lag(y.tr.num, h+1))
        
        y.fc.trend.mod <- lm(yval ~ yval.lagH + yval.lagH1, data=y.fc.trend.df)
        y.fc.newdat <- data.frame(yval.lagH = tail(y.tr.num,1), yval.lagH1 = tail(y.tr.num,2)[1])
        y.fc.trend <- predict(y.fc.trend.mod, y.fc.newdat)
      } else if (trendfc.type==4) {
        # use PCA
        
        # extract training data
        # should be of form: data frame, first column timestamp, rest: pca scores (trends!)
        pca.trend.train <- pca.trend.dat[pca.trend.dat$time <= (testsetDate-h/12+1e-8),]
        pca.trend.train <- pca.trend.train[,-1] #remove time column
        
        y.tr.num <- as.numeric(y.train.trend)
        y.fc.trend.df <- data.frame(yval = y.tr.num, 
                                    yval.lagH = lag(y.tr.num, h), 
                                    yval.lagH1 = lag(y.tr.num, h+1))
        # append all lag 2 from pca
        for (k.colname in colnames(pca.trend.train)) {
          y.fc.trend.df[paste(k.colname,1,sep="")] <- dplyr::lag(pca.trend.train[k.colname],h)
          y.fc.trend.df[paste(k.colname,2,sep="")] <- dplyr::lag(pca.trend.train[k.colname],h+1)
        }
        
        y.fc.trend.mod <- lm(yval ~ ., data=y.fc.trend.df)
        y.fc.newdat <- data.frame(yval.lagH = tail(y.tr.num,1), yval.lagH1 = tail(y.tr.num,2)[1])
        for (k.colname in colnames(pca.trend.train)) {
          y.fc.newdat[paste(k.colname,1,sep="")] <- tail(pca.trend.train[k.colname],1)
          y.fc.newdat[paste(k.colname,2,sep="")] <- tail(pca.trend.train[k.colname],2)[1,]
        }
        
        y.fc.trend <- predict(y.fc.trend.mod, y.fc.newdat)
      } else {
        y.fc.trend <- tail(y.train.trend, n=1)
      }
      
      # add trends and forecast
      y.fc <- y.fc.trend[1] + y.fc.seas + fc
    }
    
    if (useinflation) {
      y.fc <- y.fc*infl.test
      y.test <- y.test*infl.test
    }
    
    rse <- ((y.fc - y.test[1])/y.test[1])^2 
    rses[abs(time(rses) - testsetDate) < 1e-6] <- rse
    
    if (plotPredictions) {
      y.forecasts[abs(time(y.forecasts) - testsetDate) < 1e-6] <- y.fc
    }
  }
  
  if (nErr > 0) {
    print(paste0("    caught ", nErr, " error(s)"))
  }
  
  if (plotPredictions) {
    par(mfrow=c(3,1))
    plot(y.dat, col="black")
    lines(window(y.forecasts, start=start_test, end=end_test), col="red")
    
    
    plot(y.dat - y.forecasts)
    #lines(lowess(y.dat - y.forecasts, f=0.1), col="red")
    lines(computeTrend(y.dat - y.forecasts, 1, filter.size), col="red")
    
    plot(rses, type="l")
    #lines(lowess(rses, f=0.1), col="red")
    lines(computeTrend(rses, 1, filter.size), col="red")
  }
  
  # return mean relative squared error
  start_test_ind <-match(start_test, time(y.dat))
  return(list("mean"=mean(rses[-c(1:(start_test_ind-1))]), "median"=median(rses[-c(1:(start_test_ind-1))])))
}



predict_singlemonth <- function(y.dat, start_train, start_test, p = 1, q = 0, d = 0, h = 1, 
                                useinflation = FALSE, trendfc.type = 1, noisefc.type = 1, 
                                inflation.dat = NULL, plotPredictions = FALSE,
                                pca.trend.dat = NULL, pca.noise.dat = NULL) {
  
  
    testsetDate <- start_test
    end_test <- start_test
    testdates <- time(y.dat)
    start_train <- testdates[1]
  
    y.train <- y.dat
    
    
    if (useinflation) {
      infl.train <-  window(inflation.dat, start=start_train, end=testsetDate-h/12)
      y.train <- y.train/infl.train
      infl.test <- as.numeric(window(inflation.dat, start=testsetDate, end=testsetDate))
    }
    
    if (trendfc.type==-1) { # use sarima
      arima.mod <- arima(y.train, order=c(p,d,q), seasonal=list(order=c(1,0,0), period=12))
      y.fc <- predict(arima.mod,1)$pred
    } else {
      
    # detrend training data
    filter.size <- 12
    y.train.trend <- computeTrend(y.train, 1, filter.size)
    y.train.detrended <- y.train - y.train.trend
    y.train.seasyr <- computeSeas(y.train.detrended)
    y.train.rand <- removeSeas(y.train.detrended, y.train.seasyr)
    
    # fit model and forecast
    if (noisefc.type == 0) {
      fc <- 0
    } else if (noisefc.type == 1) {
      fc <- Inf
      result = tryCatch({
        fc <- forecast(Arima(y.train.rand, order=c(p,0,q)), h = h)$mean[1]
      }, error = function(e) {
        
      })
    } else if (noisefc.type == 2) {
      # use PCA predictors
      
      # build pc.df
      pc.noise.train <- pca.noise.dat[pca.noise.dat$time <= (testsetDate-h/12+1e-8),]
      
      pc.df <- data.frame(time = pc.noise.train$time)
      pc.noise.train <- pc.noise.train[,-1]
      
      # append all lag 2 from pca
      for (k.colname in colnames(pc.noise.train)) {
        pc.df[paste(k.colname,1,sep="")] <- lag(pc.noise.train[k.colname],h)
        pc.df[paste(k.colname,2,sep="")] <- lag(pc.noise.train[k.colname],h+1)
      }
      pc.df <- pc.df[,-1]
      
      ax.mod <- arimax(y.train.rand, order=c(p,d,q), xreg=pc.df)
      
      newxreg <- data.frame(time=1)
      for (k.colname in colnames(pc.noise.train)) {
        newxreg[paste(k.colname,1,sep="")] <- tail(pc.noise.train[k.colname],1)
        newxreg[paste(k.colname,2,sep="")] <- tail(pc.noise.train[k.colname],2)[1,]
      }
      newxreg <- newxreg[,-1]
      
      fc <- as.numeric(predict(ax.mod, n.ahead=h, newxreg)$pred)[h]
      
    } else {
      fc <- 0
    }
    
    # if fit or forecast unsuccessful, forecast noise with 0
    nErr <- 0
    if (fc == Inf) {
      fc <- 0
      nErr <- 1
    }
    
    # compute seasonal component and trend for forecast
    y.fc.month <- round((start_test%%1)*12)
    y.fc.seas <- y.train.seasyr[y.fc.month+1]
    
    if (trendfc.type == 1) {
      y.fc.trend <- tail(y.train.trend, n=1)
    } else if (trendfc.type == 2) {
      y.fc.trend.dat <- window(y.train, start=testsetDate-(h-1+filter.size)/12, end=testsetDate-h/12)
      y.fc.trend.df <- data.frame(xval = time(y.fc.trend.dat), yval = y.fc.trend.dat)
      colnames(y.fc.trend.df) <- c("xval", "yval")
      y.fc.trend.mod <- lm(yval ~poly(xval, 2), data=y.fc.trend.df)
      y.fc.trend <- predict(y.fc.trend.mod, data.frame(xval=start_test))
    } else if (trendfc.type==3) {
      # use last 2 values from trend (after shifting by h)
      y.tr.num <- as.numeric(y.train.trend)
      y.fc.trend.df <- data.frame(yval = y.tr.num, 
                                  yval.lagH = lag(y.tr.num, h), 
                                  yval.lagH1 = lag(y.tr.num, h+1))
      
      y.fc.trend.mod <- lm(yval ~ yval.lagH + yval.lagH1, data=y.fc.trend.df)
      y.fc.newdat <- data.frame(yval.lagH = tail(y.tr.num,1), yval.lagH1 = tail(y.tr.num,2)[1])
      y.fc.trend <- predict(y.fc.trend.mod, y.fc.newdat)
    } else if (trendfc.type==4) {
      # use PCA
      
      # extract training data
      # should be of form: data frame, first column timestamp, rest: pca scores (trends!)
      pca.trend.train <- pca.trend.dat[pca.trend.dat$time <= (testsetDate-h/12+1e-8),]
      pca.trend.train <- pca.trend.train[,-1] #remove time column
      
      y.tr.num <- as.numeric(y.train.trend)
      y.fc.trend.df <- data.frame(yval = y.tr.num, 
                                  yval.lagH = lag(y.tr.num, h), 
                                  yval.lagH1 = lag(y.tr.num, h+1))
      # append all lag 2 from pca
      for (k.colname in colnames(pca.trend.train)) {
        y.fc.trend.df[paste(k.colname,1,sep="")] <- lag(pca.trend.train[k.colname],h)
        y.fc.trend.df[paste(k.colname,2,sep="")] <- lag(pca.trend.train[k.colname],h+1)
      }
      
      y.fc.trend.mod <- lm(yval ~ ., data=y.fc.trend.df)
      y.fc.newdat <- data.frame(yval.lagH = tail(y.tr.num,1), yval.lagH1 = tail(y.tr.num,2)[1])
      for (k.colname in colnames(pca.trend.train)) {
        y.fc.newdat[paste(k.colname,1,sep="")] <- tail(pca.trend.train[k.colname],1)
        y.fc.newdat[paste(k.colname,2,sep="")] <- tail(pca.trend.train[k.colname],2)[1,]
      }
      
      y.fc.trend <- predict(y.fc.trend.mod, y.fc.newdat)
    } else {
      y.fc.trend <- tail(y.train.trend, n=1)
    }
    
    # add trends and forecast
    y.fc <- y.fc.trend[1] + y.fc.seas + fc
  }
  
  if (useinflation) {
    y.fc <- y.fc*infl.test
  }
    
  
  
  if (nErr > 0) {
    print(paste0("    caught ", nErr, " error(s)"))
  }
  
  
  return(y.fc)
}


compute_rmse_trendonly <- function(y.dat, start_train, start_test, end_test, plotPredictions = FALSE) {
  # maybe improve efficiency by only recomputing trend / seasonality for new entry in training data
  
  sqsum <- 0
  
  # dates we evaluate at
  y.alltests <- window(y.dat, start=start_test, end=end_test)
  testdates <- time(y.alltests)
  ntests <- length(testdates)
  
  if (plotPredictions) {
    y.forecasts <- y.dat
    rses <- y.dat
    rses[] <- 0
  }
  
  for (testsetDate in testdates) {
    y.train <- window(y.dat, start=start_train, end=testsetDate-1/12)
    y.test <- window(y.dat, start=testsetDate, end=testsetDate)
    
    
    # detrend training data
    filter.size <- 12
    y.train.trend <- computeTrend(y.train, 1, filter.size)
    y.train.detrended <- y.train - y.train.trend
    y.train.seasyr <- computeSeas(y.train.detrended)
    
    # no forecasting
    #y.train.rand <- removeSeas(y.train.detrended, y.train.seasyr)
    fc <- 0

    # compute seasonal component and trend for forecast
    y.fc.month <- round((time(y.test)[1]%%1)*12)
    y.fc.seas <- y.train.seasyr[y.fc.month+1]
    y.fc.trend <- tail(y.train.trend, n=1)
    
    # add trends and forecast
    y.fc <- y.fc.trend[1] + y.fc.seas + fc
    rse <- ((y.fc - y.test[1])/y.test[1])^2 
    
    if (plotPredictions) {
      y.forecasts[abs(time(y.forecasts) - testsetDate) < 1e-6] <- y.fc
      rses[abs(time(y.forecasts) - testsetDate) < 1e-6] <- rse
    }
    
    # add relative squared error
    sqsum <- sqsum + rse
    
  }
  
  if (ntests != length(testdates)) {
    print(paste0("    caught ", length(testdates)-ntests, " error(s)"))
  }
  
  if (plotPredictions) {
    par(mfrow=c(3,1))
    plot(y.dat, col="black")
    lines(window(y.forecasts, start=start_test, end=end_test), col="red")
    
    
    plot(y.dat - y.forecasts)
    #lines(lowess(y.dat - y.forecasts, f=0.1), col="red")
    lines(computeTrend(y.dat - y.forecasts, 1, filter.size), col="red")
    
    plot(rses, type="l")
    #lines(lowess(rses, f=0.1), col="red")
    lines(computeTrend(rses, 1, filter.size), col="red")
  }
  
  # return mean relative squared error
  return(sqsum/ntests)
}


compute_rmse_linregTrend <- function(y.dat, start_train, start_test, end_test, p, q, plotPredictions = FALSE) {
  # maybe improve efficiency by only recomputing trend / seasonality for new entry in training data
  
  sqsum <- 0
  
  # dates we evaluate at
  y.alltests <- window(y.dat, start=start_test, end=end_test)
  testdates <- time(y.alltests)
  ntests <- length(testdates)
  
  if (plotPredictions) {
    y.forecasts <- y.dat
    rses <- y.dat
    rses[] <- 0
  }
  
  for (testsetDate in testdates) {
    y.train <- window(y.dat, start=start_train, end=testsetDate-1/12)
    y.test <- window(y.dat, start=testsetDate, end=testsetDate)
    
    
    # detrend training data
    filter.size <- 12
    y.train.trend <- computeTrend(y.train, 1, filter.size)
    y.train.detrended <- y.train - y.train.trend
    y.train.seasyr <- computeSeas(y.train.detrended)
    y.train.rand <- removeSeas(y.train.detrended, y.train.seasyr)
    
    # fit model and forecast
    fc <- Inf
    result = tryCatch({
      fc <- forecast(Arima(y.train.rand, order=c(p,0,q)), h = 1)$mean[1]
    }, error = function(e) {
      
    })
    
    # if fit or forecast unsuccessful, forecast noise with 0
    if (fc == Inf) {
      fc <- 0
    }
    # compute seasonal component and trend for forecast
    y.fc.month <- round((time(y.test)[1]%%1)*12)
    y.fc.seas <- y.train.seasyr[y.fc.month+1]
    
    
    
    # compute trend via linreg instead
    #y.fc.trend <- tail(y.train.trend, n=1)
    
    y.fc.trend.dat <- window(y.dat, start=testsetDate-filter.size/12, end=testsetDate-1/12)
    y.fc.trend.df <- data.frame(xval = time(y.fc.trend.dat), yval = y.fc.trend.dat)
    y.fc.trend.mod <- lm(yval ~poly(xval, 2), data=y.fc.trend.df)
    y.fc.trend <- predict(y.fc.trend.mod, data.frame(xval=time(y.test)[1]))
    
    
    
    # add trends and forecast
    y.fc <- y.fc.trend[1] + y.fc.seas + fc
    
    rse <- ((y.fc - y.test[1])/y.test[1])^2 
    
    if (plotPredictions) {
      y.forecasts[abs(time(y.forecasts) - testsetDate) < 1e-6] <- y.fc
      rses[abs(time(y.forecasts) - testsetDate) < 1e-6] <- rse
    }
    
    # add relative squared error
    sqsum <- sqsum + rse
  }
  
  if (ntests != length(testdates)) {
    print(paste0("    caught ", length(testdates)-ntests, " error(s)"))
  }
  
  if (plotPredictions) {
    par(mfrow=c(3,1))
    plot(y.dat, col="black")
    lines(window(y.forecasts, start=start_test, end=end_test), col="red")
    
    
    plot(y.dat - y.forecasts)
    #lines(lowess(y.dat - y.forecasts, f=0.1), col="red")
    lines(computeTrend(y.dat - y.forecasts, 1, filter.size), col="red")
    
    plot(rses, type="l")
    #lines(lowess(rses, f=0.1), col="red")
    lines(computeTrend(rses, 1, filter.size), col="red")
  }
  
  # return mean relative squared error
  return(sqsum/ntests)
}


compute_rmse_multvar <- function(y.dat, start_train, start_test, end_test, p, q, h = 1, method=1, plotPredictions = FALSE) {
  # maybe improve efficiency by only recomputing trend / seasonality for new entry in training data
  
  sqsum <- rep(0,ncol(y.dat)-1)
  
  y.dates <- y.dat[,1]
  
  # indices we evaluate at
  start_train_ind <- which(abs(y.dates - start_train) < 1e-7)
  start_test_ind <-which(abs(y.dates - start_test) < 1e-7)
  end_test_ind <- which(abs(y.dates - end_test) < 1e-7)
  
  ntests <- end_test_ind - start_test_ind + 1
  
  nUnsuccess <- 0
  
  # initialize model
  if (method == 2) {
    mod.mar <- define.model(kvar=ncol(y.dat), ar=p, ma=q, rem.var=1)
    arp <- mod.mar$ar.pattern
    map <- mod.mar$ma.pattern
  }
  
  rses <- y.dat
  rses[,-1] <- 0
  
  if (plotPredictions) {
    y.forecasts <- y.dat
  }
  
  for (testsetInd in start_test_ind:end_test_ind) {
    y.train <- y.dat[start_train_ind:(testsetInd-h),]
    y.test <- y.dat[testsetInd:testsetInd,]
    
    
    # detrend training data
    filter.size <- 12
    y.train.trend <- computeTrend_df(y.train, 1, filter.size)
    y.train.detrended <- y.train
    y.train.detrended[,-1] <- y.train[,-1] - y.train.trend[,-1]
    y.train.seasyr <- computeSeas_df(y.train.detrended, filter.size)
    y.train.rand <- removeSeas_df(y.train.detrended, y.train.seasyr)
    
    # genertate row for testing instance
    y.fc.date <- y.dates[testsetInd]
    y.fc.dates <- y.dates[(testsetInd-(h-1)):testsetInd]
    y.fc.rows <- matrix(0, h, ncol(y.dat))
    y.fc.rows[1,] <- y.fc.dates
    
    if (method == 1) {
      mod.varma <- VARMA(y.train.rand[-(1:filter.size-1),-1], p = p, q = q, include.mean = T)
      fc <- VARMApred(mod.varma, h=h)$pred
    } else if (method == 2) {
      # fit model and forecast
      #fc <- Inf
      #result = tryCatch({
      # sink("NUL")
      #res.mar <- marima(ts(y.train.rand), ar.pattern=arp, ma.pattern=map, penalty=1)
      # mar.fc <- arma.forecast(series=ts(rbind(as.matrix(y.train.rand), y.fc.row)), marima=res.mar,
      #                         nstart=testsetInd-1, nstep=1, check=FALSE)
      # fc <- mar.fc$forecasts[-1, testsetInd]
      
      res.mar <- marima(as.matrix(y.train.rand), ar.pattern=arp, ma.pattern=map, penalty=1, Check=FALSE)
      # fcinput <- rbind(as.matrix(y.train.rand[(nrow(y.train.rand-max(p))):(nrow(y.train.rand)),]), y.fc.row)
      # mar.fc <- arma.forecast(series=fcinput, marima=res.mar,
      #                         nstart=max(p), nstep=1, check=TRUE)
      
      fcinput <- rbind(as.matrix(y.train.rand), y.fc.rows)
      
      mar.fc <- arma.forecast(series=fcinput, marima=res.mar,
                              nstart=testsetInd-1, nstep=1, check=FALSE)
      fc <- mar.fc$forecasts[-1, testsetInd]
      
      
      # sink()
      
      
      if (sum(abs(fc)) > 1e6) {
        print("large fc value!")
      }
      
      fc[abs(fc) > 1e5] <- 0
      #}, error = function(e) {
      
      #})
      
      
    }
    
      
  
    # if fit or forecast unsuccessful, forecast noise with 0
    #if (length(fc) == 1) {
    #  fc <- numeric(ncol(y.dat))
    #  nUnsuccess <- nUnsuccess+1
    #} else {
      # compute seasonal component and trend for forecast
      y.fc.month <- round((y.fc.date%%1)*12)
      y.fc.seas <- y.train.seasyr[y.fc.month+1,]
      y.fc.trend <- y.train.trend[nrow(y.train.trend),-1]
      
      # add trends and forecast
      y.fc <- y.fc.trend
      y.fc <- y.fc + y.fc.seas + fc
      rse <- ((y.fc - y.test[-1])/y.test[-1])^2
      rses[testsetInd,-1] <- rse
      
      if (plotPredictions) {
        y.forecasts[testsetInd,-1] <- y.fc
        
      }
    #}
  }
  
  if (nUnsuccess > 0) {
    print(paste0("    caught ", nUnsuccess, " error(s)"))
  }
  
  if (plotPredictions) {
    par(mfrow=c(ncol(y.dat)-1,1))
    for (k in 2:ncol(y.dat)) {
      plot(rses[,1], rses[,k])
    }
    
    for (k in 2:ncol(y.dat)) {
      par(mfrow=c(3,1))
      plot(y.dat[,1], y.dat[,k], col="black")
      lines(y.dat[start_test_ind:nrow(y.dat),1], y.forecasts[start_test_ind:nrow(y.dat),k], col="red")

      plot(y.dat[,1], y.dat[,k] - y.forecasts[,k])
      lines(computeTrend(ts(y.dat[start_test_ind:nrow(y.dat),k] - y.forecasts[start_test_ind:nrow(y.dat),k], frequency = 12, start=y.dat[start_test_ind,1], end=y.dat[nrow(y.dat),1]), 1, filter.size), col="red")

      plot(rses[,1], rses[,k], type="l")
      lines(computeTrend(ts(rses[start_test_ind:nrow(y.dat),k], frequency=12, start=rses[start_test_ind,1], end=rses[nrow(rses),1]), 1, filter.size), col="red")
    }
  }
  
  # return mean relative squared error
  return(list("mean"=colMeans(rses[-c(1:(start_test_ind-1)),-1]), "median"=(rses[-c(1:(start_test_ind-1)),-1] %>% summarise(across(where(is.numeric), median)))))
}

geom.mean <-function(x){
  return(exp(mean(log(x))))
}

quotients <-function(x){
  return(exp(diff(log(x))))
}

computeTrend <- function(y.dat, filter.type, filter.size) {
  if (filter.type == 1) {
    filter.mean <- c(rep(1, filter.size))
    filter.mean <- filter.mean / sum(filter.mean)
    
    y.trend <- stats::filter(y.dat, filter=filter.mean, method="convolution", sides=1)
  } else {
    filter.gauss.sig2 <- 1.5*filter.size
    filter.gauss <- 1/sqrt(2*pi*filter.gauss.sig2) * exp(-(0:filter.size)^2/(2*filter.gauss.sig2))
    filter.gauss <- filter.gauss/sum(filter.gauss)
    
    y.trend <- stats::filter(y.dat, filter=filter.gauss, method="convolution", sides=1)
  }
  
  switch(filter.type, 
   filt.mean={
     filter.mean <- c(rep(1, filter.size))
     filter.mean <- filter.mean / sum(filter.mean)
     
     y.trend <- stats::filter(y.dat, filter=filter.mean, method="convolution", sides=1)
   },
   filt.gauss={
     filter.gauss.sig2 <- 1.5*filter.size
     filter.gauss <- 1/sqrt(2*pi*filter.gauss.sig2) * exp(-(0:filter.size)^2/(2*filter.gauss.sig2))
     filter.gauss <- filter.gauss/sum(filter.gauss)
     
     y.trend <- stats::filter(y.dat, filter=filter.gauss, method="convolution", sides=1)  
   },
   filt.reg={
     
     print('regression filter not yet implemented')
     
        
   },
   {
     print('filter.type not implemented')
   }
  )
  
  
  return(y.trend)
}


computeTrend_df <- function(y.dat, filter.type, filter.size) {
  y.dates <- y.dat[,1]  
  y.trend <- y.dat
  
  for (i in 2:ncol(y.dat)) {
    y.ts <- ts(y.dat[,i], frequency=12, start=y.dates[1], end=y.dates[length(y.dates)])
    y.trend[,i] <- computeTrend(y.ts, filter.type, filter.size)
  }
  
  return(y.trend)
}



# better version
computeSeas <- function(y.dat){
  y.seasyr <- rep(0,12)
  for (mon in 1:12) {
    y.dat.mon <- y.dat[cycle(y.dat) == mon]
    y.seasyr[mon] <- mean(y.dat.mon, na.rm=TRUE)
  }
  
  return(y.seasyr)
}

computeSeas_df <- function(y.dat, filter.size){
  y.dates <- y.dat[,1]  
  y.seasyr <- replicate(ncol(y.dat)-1,numeric(12))
  
  for (i in 2:ncol(y.dat)) {
    y.ts <- ts(y.dat[-(1:(filter.size-1)),i], frequency=12, start=y.dates[filter.size], end=y.dates[length(y.dates)])
    y.seasyr[,i-1] <- computeSeas(y.ts)
  }
  
  return(y.seasyr)
}


# better version
removeSeas <- function(y.dat, y.seasyr) {
  y.rand <- y.dat
  
  for (mon in 1:12) {
    y.rand[cycle(y.rand) == mon] <- y.dat[cycle(y.dat) == mon] - y.seasyr[mon]
  }
  
  return(y.rand)
}

removeSeas_df <- function(y.dat, y.seasyr) {
  y.dates <- y.dat[,1]
  y.rand <- y.dat
  
  for (i in 2:ncol(y.dat)) {
    y.ts <- ts(y.dat[,i], frequency=12, start=y.dates[1], end=y.dates[length(y.dates)])
    y.rand[,i] <- removeSeas(y.ts, y.seasyr[,i-1])
  }
  
  return(y.rand)
}


removeTrendsAndSeas <- function(y.dat, filter.size) {
  y.trend <- computeTrend(y.dat, 1, filter.size)
  y.detrended <- y.dat - y.trend
  y.seasyr <- computeSeas(y.detrended)
  y.rand <- removeSeas(y.detrended, y.seasyr)
  
  return(y.rand)
}

removeTrendsAndSeas_df <- function(y.dat, filter.size) {
  y.dates <- y.dat[,1]  
  y.rand <- y.dat
  for (i in 2:ncol(y.dat)) {
    
    y.ts <- ts(y.dat[,i], frequency=12, start=y.dates[1], end=y.dates[length(y.dates)])
    
    y.trend <- computeTrend(y.ts, 1, filter.size)
    y.detrended <- y.ts - y.trend
    y.seasyr <- computeSeas(y.detrended)
    y.rand[,i] <- removeSeas(y.detrended, y.seasyr)
  }

  return(y.rand)
}


compute_noise_pca_df <- function(y.dat, ctrynames, filter.size) {
  y.temp <- data.frame(y.dat[ctrynames])
  y.temp <- cbind(y.dat$TIME_PERIOD, y.temp)
  res <- removeTrendsAndSeas_df(y.temp,filter.size)
  
  for (k.ctry in ctrynames) {
    res[,k.ctry] <- as.numeric(res[,k.ctry])
  }
  res <- res[-(1:11),]
  return(prcomp(res[,-1], scale. = TRUE))
}




# find arima(x) parameters p,d,q for each country using cross validation (evaluated in parallel)
findpdq <- function(y.eiint.dat, which.countries, infl.dat, infl_start_year, trendfc = 1, noisefc = 1, maxp=5, maxd=0, maxq=5, pca.df = NULL, pca.noise.df = NULL) {
  start_train <- min(y.eiint.dat$TIME_PERIOD)
  start_test <- 2009+1/12
  
  sysdate <- Sys.Date()
  sysmonth <- as.numeric(format(sysdate, "%m"))
  sysyr <- as.numeric(format(sysdate, "%Y"))
  end_test <- sysyr + (sysmonth-1-2)/12 # -1 for counting differently
  # -2 for last available month
  
  rseType <- "median"
  h <- 2
  useinflation <- TRUE
  
  
  
  cores <- parallel::detectCores()-1
  cl <- makeSOCKcluster(cores)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min=1, max=length(which.countries), style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  tryCatch({
    result <- foreach(k = 1:length(which.countries), .options.snow=opts)%dopar%{
      source(".\\importexport_functions.R")
      
      k.ctry <- which.countries[k]
      inflation_ts <- ts(rbind(infl.dat[k.ctry], infl.dat[nrow(infl.dat), k.ctry]), 
                         frequency = 12, start = c(infl_start_year,1))
      y.ctry.ts <- ts(y.eiint.dat[k.ctry], frequency = 12, start=y.eiint.dat$TIME_PERIOD[1], end = tail(y.eiint.dat$TIME_PERIOD,1))
      
      minErr <- Inf
      bestmod <- c(-1,0,-1)
      for (p in 0:maxp) {
        for (d in 0:maxd) {
          for (q in 0:maxq) {
            curErr <- compute_rmse(y.ctry.ts, start_train, start_test, end_test, h = h, 
                                   useinflation = useinflation, inflation.dat = inflation_ts, 
                                   p=p, d=d, q=q,  trendfc.type = trendfc, noisefc.type = noisefc, plotPredictions = FALSE,
                                   pca.trend.dat=pca.df, pca.noise.dat = pca.noise.df)[[rseType]]
            
            if (curErr < minErr) {
              bestmod <- c(p,d,q)
              minErr <- curErr
            }
          }
        }
      }
      bestmod
    }
  })
  close(pb)
  stopCluster(cl)
  return(result)
}



# compute predictions of ARIMA(X) models, given parameters p, d, q.
predictpdq <- function(y.eiint.dat, which.countries, infl.dat, infl_start_year, trendfc = 1, noisefc = 1, models, pca.df = NULL, pca.noise.df = NULL) {
  start_train <- min(y.eiint.dat$TIME_PERIOD)
  start_test <- 2009+1/12
  
  sysdate <- Sys.Date()
  sysmonth <- as.numeric(format(sysdate, "%m"))
  sysyr <- as.numeric(format(sysdate, "%Y"))
  
  start_train <- min(y.eiint.dat$TIME_PERIOD)
  sysdate <- Sys.Date()
  sysmonth <- as.numeric(format(sysdate, "%m"))
  sysyr <- as.numeric(format(sysdate, "%Y"))
  start_test <- sysyr + (sysmonth-1)/12 # current month;  -1 for counting differently
  
  h <- 2
  useinflation <- TRUE
  
  
  # parallel prediction across countries
  result <- list()
  for (k in 1:length(which.countries)) {
    k
    
    k.ctry <- which.countries[k]
    inflation_ts <- ts(infl.dat[k.ctry], frequency = 12, start = c(infl_start_year,1))
    y.ctry.ts <- ts(y.eiint.dat[k.ctry], frequency = 12, start=y.eiint.dat$TIME_PERIOD[1], end = tail(y.eiint.dat$TIME_PERIOD,1))
    mod <- models[[k]] # corr?
    
    prediction <- predict_singlemonth(y.ctry.ts, start_train, start_test, p = mod[1], q = mod[3], d = mod[2], h = h, 
                                      useinflation = useinflation, trendfc.type = trendfc, noisefc.type = noisefc, 
                                      inflation.dat = inflation_ts, plotPredictions = FALSE,
                                      pca.trend.dat = pca.df, pca.noise.dat = pca.noise.df)
    result[k] <- prediction
  }
  return(result)
}








# write output for json file
writeOutput <- function(preds, ctrynamesHR) {
  predictions = list()
  output = ""
  for(k in 1:length(ctrynamesHR)){
    k.ctry <- as.character(ctrynamesHR[k])
    predictions[k.ctry] = preds[[k]]
    output=paste(output,"    \"",k.ctry,"\": ", round(as.numeric(predictions[k]), 1) ,sep="")
    if(k.ctry!="SK"){
      output=paste(output,",\n",sep="")
    }
  }
  cat(output)
}