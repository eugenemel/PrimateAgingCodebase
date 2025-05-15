# Support Functions

#' Create a life table from survival data
#'
#' @param X A data frame containing lifespan and event status columns
#' @param interval The time interval for age grouping (default: 2)
#' @return A life table data frame
lifetable <- function(X, interval = 2) {
  X <- X %>% arrange(lifespan)
  age_breaks <- seq(from = 0, to = ceiling(max(X$lifespan + interval)), interval)
  X$age_groups <- cut(X$lifespan, age_breaks, include.lowest = TRUE, right = FALSE)
  X$int.start <- age_breaks[findInterval(X$lifespan, age_breaks)]
  X$int.end <- age_breaks[findInterval(X$lifespan, age_breaks, rightmost.closed = TRUE) + 1]
  
  
   # Helper function for cumsum with NA handling
  cumsum.na <- function(x, ...) {  x[is.na(x)] <- 0; cumsum(x)  }
  
  res <- X %>% group_by(int.start, int.end, age_groups) %>% 
    summarise(lost = sum(event == 0), events = sum(event == 1))
  
  res$int.length <- res$int.end - res$int.start
  res$int.midpoint <- with(res, int.start + (int.end - int.start) / 2)
  res$enter <- rev(cumsum(rev(res$lost + res$events)))
  res$rate <- with(res, events / (enter - events / 2 - lost / 2) / int.length)  
  
 
  res$surv <- with(res, exp(-cumsum.na(rate * int.length)))
  res$sd <- sqrt(res$rate * (1 - res$rate) / res$enter)  # binomial variance n*q*(1-q)
  res$log_rate <- log(res$rate)
  res$time <- res$int.start + (res$int.end - res$int.start) / 2
  res$sexMature = min(X$sexMature); #used for supplementary material export
  
  return(res)
}

#' Fit a simple linear model
#'
#' @param df Data frame containing log_rate, time and optionally sex columns
#' @return A linear model object
lmmodel <- function(df) { 
  if (length(unique(df$sex)) > 1) {
    model <- lm(log_rate ~ sex + time + sex:time, df)
  } else { 
    model <- lm(log_rate ~ time, df)
  } 
  return(model)
}

#' Fit a robust linear model
#'
#' @param df Data frame containing log_rate, time and optionally sex columns
#' @return A robust linear model object
rlmmodel <- function(df) { 
  tmpdf <- df
  tmpdf$sex <- paste(df$sex)
  
  if ("Both" %in% df$sex) {
    model <- MASS::rlm(log_rate ~ time, df)
  } else { 
    model <- MASS::rlm(log_rate ~ sex + time + sex:time, tmpdf)
  } 
  return(model)
}

#' Add linear model predictions to data frame
#' @param df Input data frame
#' @return Data frame with lm_fitted column added
lmpredict <- function(df) {
  model <- lmmodel(df)
  df$lm_fitted <- model$fitted
  return(df)
}

#' Get tidy output from linear model
#' @param df Input data frame for model
#' @return Tidy summary of model
lmbroom <- function(df) {
  model <- lmmodel(df)
  return(broom::tidy(model))
}

#' Add robust linear model predictions to data frame
#' @param df Input data frame
#' @return Data frame with rlm_fitted column added
rlmpredict <- function(df) {
  model <- rlmmodel(df)
  df$rlm_fitted <- model$fitted
  return(df)
}

#' Get tidy output from robust linear model
#' @param df Input data frame for model
#' @return Tidy summary of model
rlmbroom <- function(df) {
  model <- rlmmodel(df)
  return(broom::tidy(model))
}

#' Gompertzian quantile function
#' 
#' @param p Probability
#' @param shape Shape parameter
#' @param rate Rate parameter
#' @return Quantile value
qgomp <- function(p, shape, rate) {
  1 / shape * log1p(-log1p(-p) * shape / rate)
}
