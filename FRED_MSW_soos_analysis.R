rm(list=ls())
require(mboost)
require(R.matlab)
require(forecast)
require(sandwich)
data<-readMat('FRED.mat')

fileprefix='FRED_MSW'
source("helper_functions.R")

##############
## Settings ##
##############

save_directory="FREDTests/"    # Add / to end or leave empty!

# List of the series in the data file to test. E.g. 1:10, c(1,2,5,7) or just a number
lock_file="seriesfile_lock.tmp"
ids_file="unstarted_series.txt"

#Training data startyear
train_startyear=1959   # Def=1959
#Test data startyear for the first model
test_startyear=1999    # Def=1999
#How many models are estimated
test_rounds=17        # Def 17/1
#Stepping (in observations) between models.
test_blocksize=12  # Def=12/204
#Minimum number of valid observations required for estimation
min_obs_required=50    # Def=50

# Do boost-OLS?
do_bols=0

#Set to 1 to store the in-sample fits of all models (otherwise only first ones stored)
store_all_insamples=0

#Set to 1 to store g_i function graphs of all models (otherwise only first ones stored)
store_all_partials=0

# Which data series to use (d.logdiff is the transformed variable)
# MSW version does the transforms within this script.
yt_data=data$d.values

# Boosting parameters
predlags=12            # Def=12
araic_maxlags=12       # Def=12
boost_maxlags=12       # Def=12 HOX!!!! boost_use_lags setting below !!!!

hmax=12                # Def=12
hrange=c(1:12)         # Def=1:12
deg_freedom=4          # Def=4
cv_folds=5             # Def=5
boost_use_lags='FULL'   # "FULL" for all 12 or "BIC" for same as with AR(BIC)
maxlearners=300        # Def 300
bols_maxlearners=100   # Def 300
spline_bctrl = boost_control(mstop = maxlearners)
bols_bctrl = boost_control(mstop = bols_maxlearners)

# Preparations..
if(test_rounds>1) {
  esttype='recursive'
} else{
  esttype='single'
}

deptype='stationary'

fileid=sprintf("%s_%s_%s_%d_%d_%d",fileprefix,esttype,deptype,train_startyear,test_startyear,(test_rounds*test_blocksize))

#Singleton dimensions will be collapsed, so..
if(test_rounds==1) {
  store_all_insamples=0
  store_all_partials=0
}

# Which series to forecast
series<-1:128

##############
## The loop ##
##############
loop_continue<-1

ts<-series[1]
series<-series[-1]

while(loop_continue) {
  
  # Create the overall data frame
  ts_data<-data_to_frame_multistep(yt_data[,ts],predlags,hmax,data$d.year,data$d.month, data$d.transform[,ts])
  
  # Find the initial training and test data start and end indexes
  #ts_len=nrow(data$d.values)
  train_startidx=min(which(ts_data$Year %in% c(train_startyear)))
  train_endidx=max(which(ts_data$Year %in% c(test_startyear-1)))
  test_startidx=min(which(ts_data$Year %in% c(test_startyear)))
  
  #m0: NULL const, alternative boosting
  #m1: NULL AR(p), alternative boosting
  #m2: NULL AR(p), alternative AR(p)+boosting

  #Space for storing the learner counts
  bspline_learners<-matrix(data=NaN,ncol=nrow(ts_data),nrow=hmax)
  
  #Const prediction for each horizon
  const_pred<-matrix(data=NaN,ncol=nrow(ts_data),nrow=hmax)
  
  # AR(BIC) prediction and lags
  bic_pred<-matrix(data=NaN,ncol=nrow(ts_data),nrow=hmax)
  bic_lags<-matrix(data=NaN,ncol=nrow(ts_data),nrow=hmax)
  
  # Two-stage boosting
  tsboost_pred<-matrix(data=NaN,ncol=nrow(ts_data),nrow=hmax)
  tsboost_learners<-matrix(data=NaN,ncol=nrow(ts_data),nrow=hmax)
  
  # Two-stage boosting with lag restriction
  tsboost_res_pred<-matrix(data=NaN,ncol=nrow(ts_data),nrow=hmax)
  tsboost_res_learners<-matrix(data=NaN,ncol=nrow(ts_data),nrow=hmax)
  
  # No extrapolation prediction
  bspline_noextra_pred<-matrix(data=NaN,ncol=nrow(ts_data),nrow=hmax)
  bspline_noextra_mask<-matrix(data=0,ncol=nrow(ts_data),nrow=hmax)
  bspline_noextra_replaced<-mat.or.vec(hmax,1)
  tsboost_noextra_pred<-matrix(data=NaN,ncol=nrow(ts_data),nrow=hmax)
  tsboost_noextra_mask<-matrix(data=0,ncol=nrow(ts_data),nrow=hmax)
  tsboost_noextra_replaced<-mat.or.vec(hmax,1)
  
  # In-sample fitted values
  if(store_all_insamples==1) {
    ar_insample_fits<-array(data=NaN,dim=c(test_rounds,hmax,nrow(ts_data)))
    bic_insample_fits<-array(data=NaN,dim=c(test_rounds,hmax,nrow(ts_data)))
    bspline_insample_fits<-array(data=NaN,dim=c(test_rounds,hmax,nrow(ts_data)))
    tsboost_insample_fits<-array(data=NaN,dim=c(test_rounds,hmax,nrow(ts_data)))
  } else {
    ar_insample_fits<-array(data=NaN,dim=c(hmax,nrow(ts_data)))
    bic_insample_fits<-matrix(data=NaN,ncol=nrow(ts_data),nrow=hmax)
    bspline_insample_fits<-array(data=NaN,dim=c(hmax,nrow(ts_data)))
    tsboost_insample_fits<-array(data=NaN,dim=c(hmax,nrow(ts_data)))
  }
  
  # Pointers for estimation rounds
  train_data_start=mat.or.vec(test_rounds,1)
  train_data_end=mat.or.vec(test_rounds,1)
  test_data_start=mat.or.vec(test_rounds,1)
  test_data_end=mat.or.vec(test_rounds,1)
  
  cat ('\nStarting time series',paste0(data$d.varname[ts]),' (',ts,') ',train_startyear,' (',maxlearners,')\n')
  pb <- txtProgressBar(min = 0, max = test_rounds*length(hrange), style = 3)
  pbDone <- 0
  setTxtProgressBar(pb,pbDone)
  
  # Ensure reproducibility
  set.seed(12345)
  
  for (tround in 0:(test_rounds-1)) {
    
    # This points to the current y_t
    yt_idx=test_startidx+tround*test_blocksize 
    
    for (h in hrange) {
      # Name of the variable that we are trying to predict
      predicted_name<-paste0('Ytp',h)
      
      # Calculate indexes where training and test data starts and ends
      train_start=train_startidx
      train_end=yt_idx-h   # We cannot use more data in estimation than up to yt_idx
      test_start=yt_idx
      test_end=yt_idx+test_blocksize-1
      
      # Remove NaN observations from the beginning (if any)
      maxlagname=paste0('Ytm',(predlags-1))
      while(is.na(ts_data[train_start,maxlagname]) && train_start<train_end) {
        train_start=train_start+1
      }
      
      if((train_end-train_start) < min_obs_required) {
        cat("Not enough observations to estimate at h=",h)
        break
      }
      
      # Store the pointers for this round
      train_data_start[tround+1]=train_startidx
      train_data_end[tround+1]=yt_idx-1
      test_data_start[tround+1]=yt_idx
      test_data_end[tround+1]=yt_idx+test_blocksize-1
      
      #Mean, min and max values
      mu=mean(ts_data[train_start:train_end,predicted_name])
      maxval=max(max(ts_data[train_start:train_end,"Ytm11"]),max(ts_data[train_start:train_end,"Ytm0"]))
      minval=min(min(ts_data[train_start:train_end,"Ytm11"]),min(ts_data[train_start:train_end,"Ytm0"]))
      
      # Const model
      const_pred[h,test_start:test_end]<-mu
      
      # Fit AR(AIC)
      llist <- empirical_linfit_aic_h(ts_data[train_start:train_end,], araic_maxlags,h)
      lmodel <- llist$model
      ar_pred_name=paste0('ar_pred_Ytp',h)
      ar_lags_name=paste0('ar_lags_Ytp',h)
      ts_data[test_start:test_end,ar_pred_name] <- predict(lmodel, ts_data[test_start:test_end,])
      ts_data[test_start:test_end,ar_lags_name] <- llist$lags
      AIC_lags <- llist$lags
      
      # Fit AR(BIC)
      llist <- empirical_linfit_bic_h(ts_data[train_start:train_end,], araic_maxlags,h)
      biclmodel <- llist$model
      bic_pred[h,test_start:test_end] <- predict(biclmodel, ts_data[test_start:test_end,])
      bic_lags[h,test_start:test_end] <- llist$lags
      BIC_lags <- llist$lags
      
      if(boost_use_lags=="FULL") {
        boostmodel_lags=boost_maxlags
      } else {
        boostmodel_lags=BIC_lags 
        if(boostmodel_lags==0) {boostmodel_lags=1}
      }
      
      # Boost spline
      bbsmodel <- empirical_boost_spline_h(ts_data[train_start:train_end,], boostmodel_lags, spline_bctrl, deg_freedom, h)
      cvr <- cross_validate_2(bbsmodel, cv_folds)
      learners <- mstop(cvr)
      #learners <- mstop(AIC(bbsmodel))
      bspline_pred_name=paste0('bspline_pred_Ytp',h)
      bspline_learners[h,test_start:test_end]<-learners
      ts_data[test_start:test_end,bspline_pred_name]<- predict(bbsmodel[learners],ts_data[test_start:test_end,])
      
      # Boost spline to the AR(BIC) residuals
      ts_data[train_start:train_end,"BIC_residual"]<-ts_data[train_start:train_end,predicted_name]-biclmodel$fitted.values
      tsmodel <- empirical_boost_residuals(ts_data[train_start:train_end,], boostmodel_lags, spline_bctrl, deg_freedom)
      cvr <- cross_validate_2(tsmodel, cv_folds)
      learners <- mstop(cvr)
      #learners <- mstop(AIC(tsmodel))
      tsboost_learners[h,test_start:test_end]<-learners
      tsboost_pred[h,test_start:test_end]<- bic_pred[h,test_start:test_end]+predict(tsmodel[learners],ts_data[test_start:test_end,])
      
      # Boost OLS
      if(do_bols) {
        bolsmodel <- empirical_boost_ols_h(ts_data[train_start:train_end,], boostmodel_lags, bols_bctrl, h)
        cvr <- cross_validate_2(bolsmodel, cv_folds)
        blstop <- mstop(cvr)
        bols_pred_name=paste0('bols_pred_Ytp',h)
        ts_data[test_start:test_end,bols_pred_name]<- predict(bolsmodel[blstop],ts_data[test_start:test_end,])
        #bols_learners[h,(yt_idx+h):(yt_idx+h+test_blocksize-1)]=blstop
      }
      # Store in-sample values (all or just first model depending on flag)
      if(store_all_insamples) {
        ar_insample_fits[(tround+1),h,train_start:train_end]=lmodel$fitted.values
        bspline_insample_fits[(tround+1),h,train_start:train_end]=bbsmodel$fitted()
        tsboost_insample_fits[(tround+1),h,train_start:train_end]=tsmodel$fitted()
        bic_insample_fits[(tround+1),h,train_start:train_end]=biclmodel$fitted.values
      } else {
        if(tround==0) {
          ar_insample_fits[h,train_start:train_end]=lmodel$fitted.values
          bspline_insample_fits[h,train_start:train_end]=bbsmodel$fitted()
          tsboost_insample_fits[h,train_start:train_end]=tsmodel$fitted()
          bic_insample_fits[h,train_start:train_end]=biclmodel$fitted.values
        }
      }

      pbDone=pbDone+1
      setTxtProgressBar(pb,pbDone)
    }
  }
  
  train_start_str<-sprintf('01.%02d.%02d',ts_data$Month[train_start],ts_data$Year[train_start])
  test_start_str<-sprintf('01.%02d.%02d',ts_data$Month[test_startidx],ts_data$Year[test_startidx])
  
  #Make matrixes out of the prediction variables
  ar_pred=matrix(ncol=nrow(ts_data),nrow=hmax)
  ar_lags=matrix(ncol=nrow(ts_data),nrow=hmax)
  bspline_pred=matrix(ncol=nrow(ts_data),nrow=hmax)
  bols_pred=matrix(ncol=nrow(ts_data),nrow=hmax)
  true_Ytph=matrix(ncol=nrow(ts_data),nrow=hmax)
  
  for (hh in hrange) {
    ar_pred_name=paste0('ar_pred_Ytp',hh)
    ar_lags_name=paste0('ar_lags_Ytp',hh)
    ar_is_name=paste0('ar_is_Ytp',hh)
    bspline_pred_name=paste0('bspline_pred_Ytp',hh)
    bspline_is_name=paste0('bspline_is_Ytp',hh)
    bols_pred_name=paste0('bols_pred_Ytp',hh)
    Ytph_name=paste0('Ytp',hh)
    ar_pred[hh,]<-t(ts_data[,ar_pred_name])
    ar_lags[hh,]<-t(ts_data[,ar_lags_name])
    bols_pred[hh,]<-t(ts_data[,bols_pred_name])
    bspline_pred[hh,]<-t(ts_data[,bspline_pred_name])
    true_Ytph[hh,]<-ts_data[,Ytph_name]
  }
  
  filename<-sprintf('%s%s_%d(%s).mat',save_directory,fileid,ts,data$d.varname[[ts]])
  writeMat(filename,
           varname=data$d.varname[[ts]],
           vartitle=data$d.vartitle[[ts]],
           year=ts_data$Year,
           month=ts_data$Month,
           true_yt=ts_data$Ytm0,
           Xt=ts_data$Xt,
           true_Ytph=true_Ytph,
           
           const_pred=const_pred,
           ar_pred=ar_pred,
           ar_lags=ar_lags,
           ar_insample_fits=ar_insample_fits,
           
           bic_pred=bic_pred,
           bic_lags=bic_lags,
           bic_insample_fits=bic_insample_fits,
           
           bspline_pred=bspline_pred,
           bspline_learners=bspline_learners,
           bspline_insample_fits=bspline_insample_fits,
           
           bspline_noextra_pred=bspline_noextra_pred,
           bspline_noextra_mask=bspline_noextra_mask,
           bspline_noextra_replaced=bspline_noextra_replaced,
           
           bols_pred=bols_pred,
           #bols_learners=bols_learners,
           
           tsboost_learners=tsboost_learners,
           tsboost_pred=tsboost_pred,
           tsboost_insample_fits=tsboost_insample_fits,
           
           tsboost_noextra_pred=tsboost_noextra_pred,
           tsboost_noextra_mask=tsboost_noextra_mask,
           tsboost_noextra_replaced=tsboost_noextra_replaced,
           
           s_test_blocksize=test_blocksize,
           s_all_insamples=store_all_insamples,
           s_all_partials=store_all_partials,
           s_deg_freedom=deg_freedom,
           s_train_startyear=train_startyear,
           s_test_startyear=test_startyear,
           s_predlags=predlags,
           s_boost_maxlags=boost_maxlags,
           s_araic_maxlags=araic_maxlags,
           s_hmax=hmax,
           s_cv_folds=cv_folds,
           s_true_test_start=test_start_str,
           s_true_train_start=train_start_str,
           s_test_rounds=test_rounds,
           s_included_bols=do_bols,
           
           train_data_start=train_data_start,
           train_data_end=train_data_end,
           test_data_start=test_data_start,
           test_data_end=test_data_end
  )
  ################################
  # Get the next series or stop. #
  ################################
  
  if(length(series)>0) {
    ts<-series[1]
    series<-series[-1]
  } else {
    loop_continue<-0
  }
}