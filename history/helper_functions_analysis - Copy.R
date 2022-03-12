#### empirical_functions.R
#
# Apufunktioita empirical_blocksoos skriptille
#
#### Changelog
#
# 2018-03-26 TimoV: Do the stationarity transform here
# 2018-03-22 TimoV: Changelog
#

# Added saving of original Y_t
data_to_frame <- function (yt, maxp, maxh, yeardata, monthdata, Xt, transform_id) {
  
  yt_stationary<-transform_singlestep(Xt,transform_id)
  T<-length(yt_stationary)	
  framelen=T-maxp+1
  
  # Calculate the start and end indexes and create the data frame
  startidx=maxp
  endidx=T
  Year=yeardata[startidx:endidx]
  tframe<-data.frame(Year)
  y_name<-'Month'
  tframe[,y_name]=monthdata[startidx:endidx]
  
  # Y_t-x variables; 
  tstart=1
  tend=T-maxp+1
  
  for (mm in ((maxp-1):0)) {
    y_name<-paste('Ytm',mm,sep="")
    tframe[,y_name]<-yt_stationary[tstart:tend]
    if(mm==0) {
      orig_name<-paste('Xt',mm,sep="")
      tframe[,orig_name]<-Xt[tstart:tend]
    }
    tstart=tstart+1
    tend=tend+1
  }
  
  # Y_t+h variables
  for (hh in (1:maxh)) {
    y_name<-paste('Ytp',hh,sep="")
    tframe[,y_name]<-yt_stationary[tstart:tend]
    tstart=tstart+1
    tend=tend+1
  }  
  
  # Make space for the predicted variables
  pred_space=mat.or.vec(framelen,1)
  pred_space[]=NaN
  tframe$const_pred=pred_space
  
  for (ft in 1:maxh) {
    ar_yname<-paste('ar_pred_Ytp',ft,sep="")
    ar_lags_yname<-paste('ar_lags_Ytp',ft,sep="")
    ar_is_yname<-paste('ar_is_Ytp',ft,sep="")
    bspline_yname<-paste('bspline_pred_Ytp',ft,sep="")
    bspline_is_yname<-paste('bspline_is_Ytp',ft,sep="")
    bols_yname<-paste('bols_pred_Ytp',ft,sep="")
    bols_is_yname<-paste('bols_is_Ytp',ft,sep="")
    tframe[,ar_yname]<-pred_space
    tframe[,bspline_yname]<-pred_space
    tframe[,bols_yname]<-pred_space
    tframe[,ar_is_yname]<-pred_space
    tframe[,bspline_is_yname]<-pred_space
    tframe[,bols_is_yname]<-pred_space
    tframe[,ar_lags_yname]<-pred_space
  }
  
  return(tframe) 
}

# Fit linear model using "Ytph" as the dependent variable.
empirical_linfit_aic <- function(trdata, pmax) {
  
  ynam <- "Ytph"
  model<-list()
  aics<-mat.or.vec(pmax+1,1)
  
  # Do zero lags case
  fmla <- as.formula(paste(ynam, "~ 1"))
  model[[1]] <- lm(fmla, data=trdata)	 
  aics[1] <- AIC(model[[1]])
  
  if (pmax>0) {
    for (i in 2:(pmax+1)) {
      xnam <- paste("Ytm", 0:(i-2), sep="")		
      fmla <- as.formula(paste(ynam, "~", paste(xnam, collapse= "+")))
      model[[i]] <- lm(fmla, data=trdata)	 
      aics[i] <- AIC(model[[i]])
    }
  }
  nwlist=list("model"=model[[which.min(aics)]], "lags"=(which.min(aics)-1))
  
  return(nwlist)
}

empirical_linfit_aic_h <- function(trdata, pmax,h) {
  
  ynam <- paste0('Ytp',h)
  model<-list()
  aics<-mat.or.vec(pmax+1,1)
  
  # Do zero lags case
  fmla <- as.formula(paste(ynam, "~ 1"))
  model[[1]] <- lm(fmla, data=trdata)	 
  aics[1] <- AIC(model[[1]])
  
  if (pmax>0) {
    for (i in 2:(pmax+1)) {
      xnam <- paste("Ytm", 0:(i-2), sep="")		
      fmla <- as.formula(paste(ynam, "~", paste(xnam, collapse= "+")))
      model[[i]] <- lm(fmla, data=trdata)	 
      aics[i] <- AIC(model[[i]])
    }
  }
  nwlist=list("model"=model[[which.min(aics)]], "lags"=(which.min(aics)-1))
  
  return(nwlist)
}

empirical_linfit_bic_h <- function(trdata, pmax,h) {
  
  ynam <- paste0('Ytp',h)
  model<-list()
  bics<-mat.or.vec(pmax+1,1)
  
  # Do zero lags case
  fmla <- as.formula(paste(ynam, "~ 1"))
  model[[1]] <- lm(fmla, data=trdata)	 
  bics[1] <- BIC(model[[1]])
  
  if (pmax>0) {
    for (i in 2:(pmax+1)) {
      xnam <- paste("Ytm", 0:(i-2), sep="")		
      fmla <- as.formula(paste(ynam, "~", paste(xnam, collapse= "+")))
      model[[i]] <- lm(fmla, data=trdata)	 
      bics[i] <- BIC(model[[i]])
    }
  }
  nwlist=list("model"=model[[which.min(bics)]], "lags"=(which.min(bics)-1))
  
  return(nwlist)
}


empirical_linfit_fix_h <- function(trdata, lags, h) {
  
  ynam <- paste0('Ytp',h)
  
  # Do zero lags case
  if(lags==0) {
    fmla <- as.formula(paste(ynam, "~ 1"))
    lmodel <- lm(fmla, data=trdata)	 
  } else {
    xnam <- paste("Ytm", 0:(lags-1), sep="")		
    fmla <- as.formula(paste(ynam, "~", paste(xnam, collapse= "+")))
    lmodel <- lm(fmla, data=trdata)	 
  }
  
  nwlist=list("model"=lmodel, "lags"=lags)
  
  return(nwlist)
}

empirical_boost_spline <- function(trdata, maxlags, bctrl, df) {
  
  bbsxnam <- paste0("bbs(Ytm", 0:(maxlags-1), ", df=",df,")")
  ynam <- "Ytph"
  
  fmla <- as.formula(paste(ynam, "~", paste(bbsxnam, collapse="+")))
  
  model<-gamboost(fmla, data=trdata, family=Gaussian(), control = bctrl)
  
  return(model)
}

empirical_boost_spline_h <- function(trdata, maxlags, bctrl, df, h, loss_family=Gaussian()) {
  
  bbsxnam <- paste0("bbs(Ytm", 0:(maxlags-1), ", df=",df,")")
  ynam <- paste0('Ytp',h)
  
  fmla <- as.formula(paste(ynam, "~", paste(bbsxnam, collapse="+")))
  
  model<-gamboost(fmla, data=trdata, family=loss_family, control = bctrl)
  
  return(model)
}


empirical_boost_residuals <- function(trdata, maxlags, bctrl, df, loss_family=Gaussian()) {
  
  bbsxnam <- paste0("bbs(Ytm", 0:(maxlags-1), ", df=",df,")")
  ynam <- "BIC_residual"
  
  fmla <- as.formula(paste(ynam, "~", paste(bbsxnam, collapse="+")))
  
  model<-gamboost(fmla, data=trdata, family=loss_family, control = bctrl)
  
  return(model)
}

empirical_boost_residuals_bols <- function(trdata, maxlags, bctrl) {
  
  bbsxnam <- paste0("bols(Ytm", 0:(maxlags-1), ", df=",df,")")
  ynam <- "BIC_residual"
  
  fmla <- as.formula(paste(ynam, "~", paste(bbsxnam, collapse="+")))
  
  model<-gamboost(fmla, data=trdata, family=Gaussian(), control = bctrl)
  
  return(model)
}


empirical_boost_spline_xt <- function(trdata, maxlags, bctrl, df) {
  
  bbsxnam <- paste0("bbs(Ytm", 0:(maxlags-1), ", df=",df,")")
  bbsxtnam <- paste0("bbs(Xtm", 0:(maxlags-1), ", df=",df,")")

  ynam <- "Ytph"
  
  fmla <- as.formula(paste(ynam, "~", paste0(paste(bbsxnam, collapse="+"),"+",paste(bbsxtnam, collapse="+"))))
  
  model<-gamboost(fmla, data=trdata, family=Gaussian(), control = bctrl)
  
  return(model)
}

# Argument "lags" should include e.g. c(0,2,5)
empirical_boost_spline_lags <- function(trdata, lags, bctrl, df) {
  
  bbsxnam <- paste0("bbs(Ytm", lags, ", df=",df,")")
  ynam <- "Ytph"
  
  fmla <- as.formula(paste(ynam, "~", paste(bbsxnam, collapse="+")))
  
  model<-gamboost(fmla, data=trdata, family=Gaussian(), control = bctrl)
  
  return(model)
}

empirical_linfit <- function(trdata, pmax, lags=NULL) {
  
  ynam <- "Ytph"
  model<-list()
  bics<-mat.or.vec(pmax+1,1)
  
  if(is.null(lags)) {
    # Do zero lags case
    fmla <- as.formula(paste(ynam, "~ 1"))
    model[[1]] <- lm(fmla, data=trdata)	 
    bics[1] <- BIC(model[[1]])
    
    if (pmax>0) {
      for (i in 2:(pmax+1)) {
        xnam <- paste("Ytm", 0:(i-2), sep="")		
        fmla <- as.formula(paste(ynam, "~", paste(xnam, collapse= "+")))
        model[[i]] <- lm(fmla, data=trdata)	 
        bics[i] <- BIC(model[[i]])
      }
    }
    nwlist=list("model"=model[[which.min(bics)]], "lags"=(which.min(bics)-1))
  }
  else {
    xnam <- paste("Ytm", 0:(lags-1), sep="")		
    fmla <- as.formula(paste(ynam, "~", paste(xnam, collapse= "+")))
    model <- lm(fmla, data=trdata)	 
    nwlist=list("model"=model, "lags"=lags)
  }
  return(nwlist)
}

empirical_linfit_h <- function(trdata, pmax, lags=NULL, h) {
  
  ynam <- paste0("Ytp",h)
  model<-list()
  bics<-mat.or.vec(pmax+1,1)
  
  if(is.null(lags)) {
    # Do zero lags case
    fmla <- as.formula(paste(ynam, "~ 1"))
    model[[1]] <- lm(fmla, data=trdata)	 
    bics[1] <- BIC(model[[1]])
    
    if (pmax>0) {
      for (i in 2:(pmax+1)) {
        xnam <- paste("Ytm", 0:(i-2), sep="")		
        fmla <- as.formula(paste(ynam, "~", paste(xnam, collapse= "+")))
        model[[i]] <- lm(fmla, data=trdata)	 
        bics[i] <- BIC(model[[i]])
      }
    }
    nwlist=list("model"=model[[which.min(bics)]], "lags"=(which.min(bics)-1))
  }
  else {
    xnam <- paste("Ytm", 0:(lags-1), sep="")		
    fmla <- as.formula(paste(ynam, "~", paste(xnam, collapse= "+")))
    model <- lm(fmla, data=trdata)	 
    nwlist=list("model"=model, "lags"=lags)
  }
  return(nwlist)
}

empirical_boost_ols <- function(trdata, maxlags, bctrl) {
  
  xnam <- paste0(paste("Ytm", 0:(maxlags-1), sep=""))
  ynam <-"Ytph"
  fmla <- as.formula(paste(ynam, "~", paste(xnam, collapse="+")))
  
  model<-gamboost(fmla, data=trdata, family=Gaussian(), baselearner="bols", control = bctrl)
  
  return(model)
}

empirical_boost_ols_h <- function(trdata, maxlags, bctrl, h) {
  
  xnam <- paste0(paste("Ytm", 0:(maxlags-1), sep=""))
  ynam <-paste0('Ytp',h)
  fmla <- as.formula(paste(ynam, "~", paste(xnam, collapse="+")))
  
  model<-gamboost(fmla, data=trdata, family=Gaussian(), baselearner="bols", control = bctrl)
  
  return(model)
}

## MSW version related stuff

# Do the selected FRED-MD transform to yt
transform_singlestep <- function (yt, transform_id) {
  
  ytlen=length(yt)
  transformed=mat.or.vec(ytlen,1)
  transformed[]=NA
  
  if(transform_id==1) {
    transformed=yt
  }
  if(transform_id==2) {
    transformed[2:ytlen]=yt[2:ytlen]-yt[1:(ytlen-1)]
  }
  if(transform_id==3) {
    transformed[3:ytlen]=(yt[3:ytlen]-yt[2:(ytlen-1)])-(yt[2:(ytlen-1)]-yt[1:(ytlen-2)])
  }
  if(transform_id==4) {
    transformed=log(yt)
  }
  if(transform_id==5) {
    transformed[2:ytlen]=log(yt[2:ytlen])-log(yt[1:(ytlen-1)])
  }
  if(transform_id==6) {
    transformed[3:ytlen]=(log(yt[3:ytlen])-log(yt[2:(ytlen-1)]))-(log(yt[2:(ytlen-1)])-log(yt[1:(ytlen-2)]))
  }
  if(transform_id==7) {
    transformed[2:ytlen]=yt[2:ytlen]/yt[1:(ytlen-1)]-1
  }
  return(transformed)
}

# Form Y_t+h == y_t+h-y_t according to FRED-MD transform for a given h
transform_multistep <- function (yt, h, transform_id) {
  
  ytlen=length(yt)
  transformed=mat.or.vec(ytlen,1)
  transformed[]=NA
  
  if(transform_id==1) {
    transformed[1:(ytlen-h)]=yt[(h+1):ytlen]
  }
  if(transform_id==2) {
    transformed[1:(ytlen-h)]=yt[(h+1):ytlen]-yt[1:(ytlen-h)]
  }
  if(transform_id==3) {
    transformed[2:(ytlen-h)]=(yt[(h+2):(ytlen)]-yt[2:(ytlen-h)])-h*(yt[((2):(ytlen-h))]-yt[1:(ytlen-h-1)])
  }
  if(transform_id==4) {
    transformed[1:(ytlen-h)]=log(yt[(h+1):ytlen])
  }
  if(transform_id==5) {
    transformed[1:(ytlen-h)]=log(yt[(h+1):ytlen])-log(yt[1:(ytlen-h)])
  }
  if(transform_id==6) {
    transformed[2:(ytlen-h)]=(log(yt[(h+2):(ytlen)])-log(yt[2:(ytlen-h)]))-h*(log(yt[((2):(ytlen-h))])-log(yt[1:(ytlen-h-1)]))
  }
  if(transform_id==7) {
    transformed[1:(ytlen-h)]=yt[(h+1):ytlen]/yt[1:(ytlen-h)]-1
  }
  return(transformed)
}

# Leave yt as level or log-level depending on transform_id
transform_level <- function (yt, transform_id) {
  
  ytlen=length(yt)
  transformed=mat.or.vec(ytlen,1)
  transformed[]=NA
  
  if(transform_id==1) {
    transformed=yt
  }
  if(transform_id==2) {
    transformed=yt
  }
  if(transform_id==3) {
    transformed=yt
  }
  if(transform_id==4) {
    transformed=log(yt)
  }
  if(transform_id==5) {
    transformed=log(yt)
  }
  if(transform_id==6) {
    transformed=log(yt)
  }
  if(transform_id==7) {
    transformed=yt
  }
  return(transformed)
}


# Transform yt according to the transform_id. 
data_to_frame_multistep <- function (yt, maxp, maxh, yeardata, monthdata, transform_id) {
  
  yt_stationary<-transform_singlestep(yt,transform_id)
  T<-length(yt_stationary)	
  
  # Calculate the start and end indexes and create the data frame
  tframe<-data.frame(Year=t(yeardata))
  tframe$Month=t(monthdata)
  
  # Y_t-x variables
  for (mm in ((maxp-1):0)) {
    y_name<-paste('Ytm',mm,sep="")
    tframe[,y_name]<-lagpad(yt_stationary,mm)
    if(mm==0) {
      orig_name='Xt'
      tframe[,orig_name]=yt
    }
  }
  
  # Y_t+h variables and Y_t-h multistep variables
  for (hh in (1:maxh)) {
    Ytph<-transform_multistep(yt,hh,transform_id)
    y_name<-paste('Ytp',hh,sep="")
    tframe[,y_name]<-Ytph
    y_name<-paste('YtD',hh,sep="")
    tframe[,y_name]<-lagpad(Ytph,hh)
  }  
  
  # Make space for the predicted variables
  pred_space=mat.or.vec(T,1)
  pred_space[]=NA
  tframe$const_pred=pred_space
  
  for (ft in 1:maxh) {
    ar_yname<-paste('ar_pred_Ytp',ft,sep="")
    ar_lags_yname<-paste('ar_lags_Ytp',ft,sep="")
    ar_is_yname<-paste('ar_is_Ytp',ft,sep="")
    bspline_yname<-paste('bspline_pred_Ytp',ft,sep="")
    bspline_is_yname<-paste('bspline_is_Ytp',ft,sep="")
    bols_yname<-paste('bols_pred_Ytp',ft,sep="")
    bols_is_yname<-paste('bols_is_Ytp',ft,sep="")
    tframe[,ar_yname]<-pred_space
    tframe[,bspline_yname]<-pred_space
    tframe[,bols_yname]<-pred_space
    tframe[,ar_is_yname]<-pred_space
    tframe[,bspline_is_yname]<-pred_space
    tframe[,bols_is_yname]<-pred_space
    tframe[,ar_lags_yname]<-pred_space
  }
  
  return(tframe) 
}

# Cross validation for boosting models
cross_validate_2 <- function(model, folds) {
  
  cvm <- cv(model.weights(model), type = "kfold", B=folds)
  cvr <- cvrisk(model, folds=cvm, papply = lapply)
  return(cvr)
}

# Transform yt+h according to the transform_id. 
# Explanatory variables in (log) levels.
data_to_frame_level <- function (yt, maxp, maxh, yeardata, monthdata, transform_id) {
  
  yt_level<-transform_level(yt,transform_id)
  T<-length(yt_level)	
  framelen=T-maxp+1
  
  # Calculate the start and end indexes and create the data frame
  startidx=maxp
  endidx=T
  Year=yeardata[startidx:endidx]
  tframe<-data.frame(Year)
  y_name<-'Month'
  tframe[,y_name]=monthdata[startidx:endidx]
  
  # Y_t-x variables; 
  tstart=0
  tend=T-maxp
  
  for (mm in ((maxp-1):0)) {
    tstart=tstart+1
    tend=tend+1
    y_name<-paste('Ytm',mm,sep="")
    tframe[,y_name]<-yt_level[tstart:tend]
    if(mm==0) {
      orig_name='Xt'
      tframe[,orig_name]=yt[tstart:tend]
    }
  }
  
  # Y_t+h variables
  for (hh in (1:maxh)) {
    tstart=tstart+1
    tend=tend+1    
    y_name<-paste('Ytp',hh,sep="")
    tframe[,y_name]<-yt_level[tstart:tend]
  }  
  
  # Make space for the predicted variables
  pred_space=mat.or.vec(framelen,1)
  pred_space[]=NA
  tframe$const_pred=pred_space
  
  for (ft in 1:maxh) {
    ar_yname<-paste('ar_pred_Ytp',ft,sep="")
    ar_lags_yname<-paste('ar_lags_Ytp',ft,sep="")
    ar_is_yname<-paste('ar_is_Ytp',ft,sep="")
    bspline_yname<-paste('bspline_pred_Ytp',ft,sep="")
    bspline_is_yname<-paste('bspline_is_Ytp',ft,sep="")
    bols_yname<-paste('bols_pred_Ytp',ft,sep="")
    bols_is_yname<-paste('bols_is_Ytp',ft,sep="")
    tframe[,ar_yname]<-pred_space
    tframe[,bspline_yname]<-pred_space
    tframe[,bols_yname]<-pred_space
    tframe[,ar_is_yname]<-pred_space
    tframe[,bspline_is_yname]<-pred_space
    tframe[,bols_is_yname]<-pred_space
    tframe[,ar_lags_yname]<-pred_space
  }
  
  return(tframe) 
}

# Lag series x by k and pad with NA
lagpad <- function(x, k) {
  if (!is.vector(x)) 
    stop('x must be a vector')
  if (!is.numeric(x)) 
    stop('x must be numeric')
  if (!is.numeric(k))
    stop('k must be numeric')
  if (1 != length(k))
    stop('k must be a single number')
  c(rep(NA, k), x)[1 : length(x)] 
}

# Lead series x by k and pad with NA
leadpad <- function(x, k) {
  if (!is.vector(x)) 
    stop('x must be a vector')
  if (!is.numeric(x)) 
    stop('x must be numeric')
  if (!is.numeric(k))
    stop('k must be numeric')
  if (1 != length(k))
    stop('k must be a single number')
  c(x[(k+1):(length(x))], rep(NA, k))[1 : length(x)] 
}



# Lead series x by k and pad with NA
back_to_level = function(X, Y) {
  z = c(NA)
  for (i in 1:(length(X)-1)) {
    z = c(z, Y[i]*exp(X[i+1]))
  }
  return(z)
}

mape <- function(X, Y) {
  seq_an = !is.na(X)
  ind = 100*sum(abs((Y[seq_an] - X[seq_an])/Y[seq_an]))
  return(ind)
}

mae <- function(X, Y) {
  seq_an = !is.na(X)
  ind = mean(abs(Y[seq_an] - X[seq_an]))
  return(ind)
}

mse <- function(X, Y) {
  seq_an = !is.na(X)
  ind = mean((Y[seq_an] - X[seq_an])^2)
  return(ind)
}

rmse <- function(X, Y) {
  seq_an = !is.na(X)
  ind = sqrt(mean((Y[seq_an] - X[seq_an])^2))
  return(ind)
}