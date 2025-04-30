# This code is licensed under the MIT License. See the LICENSE file for details.


# This is the main file for running the cross-validation for the ARMA and ARMAX models.
# First, the data is imported and preprocessed.
# Then, the cross-validation is run to determine the optimal models for the ARMA and ARMAX models.
# Finally, the predictions are computed we produce output that can be directly pasted to relevant json files.




source(".\\importexport_data.R")


# run cross-validation to determine optimal models (ARMA)
# beep when done
sysdate <- Sys.Date()
sysmonth <- as.numeric(format(sysdate, "%m"))
basepath <- paste(".\\saved_variables_julius\\", sysmonth, "_", sep="")

res.exp.arma <- findpdq(y.exp.all, ctrynamesHR, infl.dat, infl_start_year, 
                        trendfc = 1, noisefc = 1)
saveRDS(res.exp.arma, paste(basepath, "exp.arma.rds", sep=""))
beep(5)

res.imp.arma <- findpdq(y.imp.all, ctrynamesHR, infl.dat, infl_start_year, 
                        trendfc = 1, noisefc = 1)
saveRDS(res.imp.arma, paste(basepath, "imp.arma.rds", sep=""))
beep(5)


res.exp.int.arma <- findpdq(y.exp.int.all, ctrynamesHR, infl.dat, infl_start_year, 
                            trendfc = 1, noisefc = 1)
saveRDS(res.exp.int.arma, paste(basepath, "exp.int.arma.rds", sep=""))
beep(5)



# run cross-validation to determine optimal models (ARMAY)
# beep when done
res.exp.armax <- findpdq(y.exp.all, ctrynamesHR, infl.dat, infl_start_year, 
                         trendfc = 4, noisefc = 2, pca.df = pca.df, pca.noise.df = pca.noise.df)
saveRDS(res.exp.armax, paste(basepath, "exp.armax.rds", sep=""))
beep(5)




res.imp.armax <- findpdq(y.imp.all, ctrynamesHR, infl.dat, infl_start_year, 
                         trendfc = 4, noisefc = 2, pca.df = pca.df, pca.noise.df = pca.noise.df, maxd=0)
saveRDS(res.imp.armax, paste(basepath, "imp.armax.rds", sep=""))
beep(5)




res.exp.int.armax <- findpdq(y.exp.int.all, ctrynamesHR, infl.dat, infl_start_year, 
                             trendfc = 4, noisefc = 2, pca.df = pca.df, pca.noise.df = pca.noise.df, maxd=0)
saveRDS(res.exp.int.armax, paste(basepath, "exp.int.armax.rds", sep=""))
beep(5)

beep(11)


# in case the workspace is gone, reload the model parameters.
res.exp.arma <- readRDS(paste(basepath, "exp.arma.rds", sep=""))
res.imp.arma <- readRDS(paste(basepath, "imp.arma.rds", sep=""))
res.exp.int.arma <- readRDS(paste(basepath, "exp.int.arma.rds", sep=""))
res.exp.armax <- readRDS(paste(basepath, "exp.armax.rds", sep=""))
res.imp.armax <- readRDS(paste(basepath, "imp.armax.rds", sep=""))
res.exp.int.armax <- readRDS(paste(basepath, "exp.int.armax.rds", sep=""))


# compute predictions
pred.arma.exp <- predictpdq(y.exp.all, ctrynamesHR, infl.dat, infl_start_year, trendfc = 1, noisefc = 1, models=res.exp.arma)
pred.arma.imp <- predictpdq(y.imp.all, ctrynamesHR, infl.dat, infl_start_year, trendfc = 1, noisefc = 1, models=res.imp.arma)
pred.arma.exp.int <- predictpdq(y.exp.int.all, ctrynamesHR, infl.dat, infl_start_year, trendfc = 1, noisefc = 1, models=res.exp.int.arma)


pred.armax.exp <- predictpdq(y.exp.all, ctrynamesHR, infl.dat, infl_start_year, trendfc = 4, noisefc = 2, models=res.exp.armax,  pca.df = pca.df, pca.noise.df = pca.noise.df)
pred.armax.imp <- predictpdq(y.imp.all, ctrynamesHR, infl.dat, infl_start_year, trendfc = 4, noisefc = 2, models=res.imp.armax,  pca.df = pca.df, pca.noise.df = pca.noise.df)
pred.armax.exp.int <- predictpdq(y.exp.int.all, ctrynamesHR, infl.dat, infl_start_year, trendfc = 4, noisefc = 2, models=res.exp.int.armax,  pca.df = pca.df, pca.noise.df = pca.noise.df)

beep(5)

# check for NA-values
sum(is.na(c(pred.arma.exp, pred.arma.imp, pred.arma.exp.int, pred.armax.exp, pred.armax.imp, pred.armax.exp.int)))
beep()



# output for json files

# folder exp-eu27
writeOutput(pred.arma.exp.int, ctrynamesHR) # entry 3
writeOutput(pred.armax.exp.int, ctrynamesHR) # entry 4

# folder exp-ext-eu27
writeOutput(pred.arma.exp, ctrynamesHR)
writeOutput(pred.armax.exp, ctrynamesHR)

# folder imp-ext-eu27
writeOutput(pred.arma.imp, ctrynamesHR)
writeOutput(pred.armax.imp, ctrynamesHR)








