# helper functions definition
library(KFAS)

calc_log_likelihood <- function(kfs){
    # calculate log likelihood
    #  NOTE: second term: modify term (Durbin and Koopman,2012)
    ll <- kfs$logLik - sum(kfs$Finf>0) * log(2*pi)/2
    return(ll)
}

calc_AIC_from_ll <- function(ll,mod){
    # calculate AIC using log likelihood & model
    #     AIC = -2 * maximum likelihood +
    #     2 * (number of unknown params + number of dynamic states)

    # number of params
    num_p <- sum(is.na(mod$Q)) + 1
    # number of states
    num_s <- length(mod$a1)
    aic <- -2*ll + 2*(num_p + num_s)
    return(aic)
}

calc_AIC <- function(kfs,mod){
    # calculate AIC
    # input :
    #  kfs : KFS object
    #  mod : model before parameters optimization
    ll <- calc_log_likelihood(kfs)
    aic <- calc_AIC_from_ll(ll,mod)
    return(aic)
}

calc_AICs <- function(kfss,mods){
  n_mods <- length(mods)
  aics <- vector()
  for(i in 1:n_mods){
    aics[i] <- calc_AIC(kfs[[i]],mods[[i]])
  }
  return(aics)
}

select_best_model <- function(kfss,mods){
  # return best model which minimize AIC
  n_mods <- length(mods)

  aics <- vector()
  for(i in 1:n_mods){
    aic <- calc_AIC(kfss[[i]],mods[[i]])
    aics[i] <- aic
  }
  best_idx <- which(aics == min(aics))
  return(list(best_index=best_idx,model=mods[[best_idx]]))
}

get_standardized_error <- function(kfs){
  # get standardized error

  obs_err <- rstandard(kfs,"pearson")
  state_err <- rstandard(kfs,"state")
  return(cbind(obs_err,state_err))
}


fit_wrapper <- function(mod,method="BFGS"){
  num_p <- sum(is.na(mod$Q)) + 1
  fit <- fitSSM(mod,numeric(num_p),method=method)
  return(fit)
}

get_models <- function(fits){
  n_fits <- length(fits)
  mods <- list()
  for(i in 1:n_fits){
    mods[[i]] <- fits[[i]]$model
  }
  return(mods)
}

get_calendar <- function(start,end){
  dates <- seq(as.Date(start),as.Date(end),by=1)
  weeks <- table(substr(dates,1,7),weekdays(dates,T))
  sun <- weeks[,"日"]
  mon <- weeks[,"月"]-sun; tue <- weeks[,"火"]-sun; wed <- weeks[,"水"]-sun
  thu <- weeks[,"木"]-sun; fry <- weeks[,"金"]-sun; sat <- weeks[,"土"]-sun
  calendar <- cbind(mon, tue, wed, thu, fry, sat)


  monthyear <- rownames(weeks)

  leapyear <- rep(F,length(monthyear))
  years <- as.numeric(substr(monthyear,1,4))
  feb_idx <- grep("-02",monthyear)
  not_leapyear_idx <- which(years%%4 != 0)

  leapyear[feb_idx] <- T
  leapyear[not_leapyear_idx] <- F

  return(list(calendar=calendar,leapyear=leapyear))
}

## メモ
# 予測：predict("confidence")
# ローカルレベル推定量：kfs$alphahat[,"level"]
# スロープ推定量：kfs$alphahat[,"slope"]
# error: rstandard
