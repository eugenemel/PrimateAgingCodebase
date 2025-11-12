# Structural Bayesian Regression of Body Weight and Lifespan
# Using brms with a custom family based on the Gompertz-Makeham survival function

# input brmsdata has columns:
#  Log2BodyWt = log2(MedianBodyWeight)
#  maxLifespan = observed maximum lifespan
#  species = species name
#  alpha = shape parameter of Gompertz-Makeham fitted to survival data
#  beta = rate parameter of Gompertz-Makeham fitted to survival data
#  p = quantile to be modeled (set to 0.99 here)
#  A = phylogenetic covariance matrix

#regresion model structure:
#  beta is function of body weight
#  alpha is function of body weight
#  maxLifespan  is function of alpha and beta where p=0.99
# regress maxLifespan as function alpha and beta using the qgomp2 function
#use chain of equations to implement custom family in brms
#ignore phylogenetic effect for now
# beta is between 0 and 1
# alpha is always negative
# maxLifespan is  always positive


library(brms)
#lambda BW ~ beta 0.9772544
#lambda BW ~ alpha 0.9838139
A <- ape::vcv.phylo(primate_tree_subset, model="lambda", value=0.98) #model = "lambda" for pagel` lambda

brmsdata=XIMP; 
brmsdata$species=row.names(XIMP)
brmsdata$logalpha = brmsdata$alpha
qmax = 0.99;  # quantile to model

qgomp2 <- function(p, beta, logalpha) { 
  beta=pmax(pmin(beta, 0.99), 0.001);
  logalpha=pmax(pmin(logalpha,-1), -20);
  1 / beta * log1p(-log1p(-p) * beta /exp(logalpha));
 }

# Define posterior prediction function
posterior_predict_qgomp <- function(i, prep, ...) {
  logalpha <- prep$dpars$logalpha[, i]
  beta <- prep$dpars$beta[, i]
  p <- 0.99  # Fixed quantile level
  
  # Expected quantile from Gompertz
  expected_quantile <- qgomp2(p, beta, logalpha)

  # Generate predictions with normal noise
  rnorm(prep$ndraws, expected_quantile, 0.1)
}

# Define posterior_epred function
posterior_epred_gompertz_quantile <- function(prep) {
  logalpha <- prep$dpars$logalpha
  beta <- prep$dpars$beta
  p <- 0.99  # Fixed quantile level
  
  # Expected quantile from Gompertz
  qgomp2(p, beta, logalpha)
}


# Define the custom family for maxLifespan
qgomp_family <- custom_family( "qgomp", 
  dpars = c("mu", "beta", "logalpha"),
  links = c("identity", "identity", "identity"),
  lb = c(0, 0.001, -20),
  ub = c(NA, 0.99, -1),
  type = "real",
  posterior_predict = posterior_predict_qgomp,
  posterior_epred = posterior_epred_gompertz_quantile
)

# Stan functions for the custom family

stan_funs <- "
  // Define qgomp2 function for Stan
  real qgomp2(real p, real beta, real logalpha) {
    return 1.0 / beta * log1p(-log1p(-p) * beta / exp(logalpha));
  }
  
  real qgomp_lpdf(real y, real mu, real beta, real logalpha) {
    return student_t_lpdf(y | 3, mu, 1); #centered 
  }
  
  //randomnumber generator for qgomp
  real qgomp_rng(real mu, real beta, real logalpha) {
    real u = uniform_rng(0, 1); 
    return qgomp2(u, beta, logalpha);
  }

  // functions for posterior predictions
  real posterior_pred_qgomp(real beta, real logalpha) {
    return qgomp2(0.99, beta, logalpha);
  }
"

# Stan variables needed for the custom family
stanvars <- stanvar(scode = stan_funs, block = "functions")

# To include phylogenetic effects, the model would be modified as follows:
bf_beta_phylo <- bf(beta ~ 1+Log2BodyWt + (1|gr(species, cov = A)))
bf_logalpha_phylo <- bf(logalpha ~ 1+Log2BodyWt + (1|gr(species, cov = A)))
bf_lifespan <- bf(maxLifespan ~ qgomp2(0.99, beta, logalpha), family = qgomp_family)

fit_gomp_phylo <- brm(
  bf_beta_phylo + bf_logalpha_phylo + bf_lifespan,
  data = brmsdata,
  data2 = list(A = A),
  stanvars = stanvars,
  chains = 4,
  cores = 4,
  iter = 4000, 
  control = list(adapt_delta = 0.95),
);


### plot fitted results
library(ggplot2)
library(tidybayes)

PLOTS=plot(conditional_effects(fit_gomp_phylo, effects = "Log2BodyWt"), 
           re_formula = NA, points = TRUE, rug = TRUE,ask=F,plot=F)
ggarrange(plotlist = PLOTS,nrow = 1)


#prediction
# Extract posterior draws
post <- as_draws_df(fit_gomp_phylo)

# To visualize the effects of alpha and beta on lifespan, we can generate
# predictions from the model over a grid of values.

if(1) { 
pdf("plots/body_weight_sem.pdf", width = 5, height = 12)
par(mfrow = c(3, 1), cex.lab = 1.2, cex.axis = 1.3, cex.main = 1.4)

#extract estimates
coefout=post %>%
  summarise(
    beta_Intercept_50 = quantile(b_beta_Intercept, 0.5),
    beta_Log2BodyWt_50 = quantile(b_beta_Log2BodyWt, 0.5),
    logalpha_Intercept_50 = quantile(b_logalpha_Intercept, 0.5),
    logalpha_Log2BodyWt_50 = quantile(b_logalpha_Log2BodyWt, 0.5),
) 

#predict upper and lower bounds of beta, logalpha, maxLifespan

predout=XIMP %>% arrange(Log2BodyWt)
predout$pred_beta = coefout$beta_Intercept_50+coefout$beta_Log2BodyWt_50*predout$Log2BodyWt;
predout$pred_alpha = coefout$logalpha_Intercept_50+ coefout$logalpha_Log2BodyWt_50*predout$Log2BodyWt;
predout$pred_maxLifespan = qgomp2(0.99,predout$pred_beta,predout$pred_alpha)

  #random sample 100 from post
  postsubset = post %>% sample_n(100)
  plot(predout$Log2BodyWt, predout$beta, pch=19, main="Effect of Body Weight on Aging Rate", xlab="Log2(Body Weight)", ylab="Predicted beta");
  lines(predout$Log2BodyWt, predout$pred_beta, col='#3D6FB6', lwd=5);
  #add points from posterior samples
  for(i in 1:nrow(postsubset)){
    samp_beta = postsubset$b_beta_Intercept[i]+postsubset$b_beta_Log2BodyWt[i]*predout$Log2BodyWt
    lines(predout$Log2BodyWt, samp_beta, col = "#3D6FB605",lwd=10);
  }

  plot(predout$Log2BodyWt, predout$alpha, pch=19, main="Effect of Body Weight on Baseline Hazard", xlab="Log2(Body Weight)", ylab="Predicted ln(alpha)");
  lines(predout$Log2BodyWt, predout$pred_alpha, col='#3D6FB6', lwd=5);
  #add points from posterior samples
  for(i in 1:nrow(postsubset)){
    samp_alpha = postsubset$b_logalpha_Intercept[i]+postsubset$b_logalpha_Log2BodyWt[i]*predout$Log2BodyWt
    lines(predout$Log2BodyWt, samp_alpha, col = "#3D6FB605",lwd=10);
  }
  
  plot(predout$Log2BodyWt, predout$maxLifespan, pch=19, main="Effect of Body Weight on Max.Lifespan", xlab="Log2(Body Weight)", ylab="Predicted Max Lifespan");
  lines(predout$Log2BodyWt, predout$pred_maxLifespan, col="#3D6FB6", lwd=5);
  #add points from posterior samples
  for(i in 1:nrow(postsubset)){
    samp_lifespan = qgomp2(0.99,
      postsubset$b_beta_Intercept[i]+postsubset$b_beta_Log2BodyWt[i]*predout$Log2BodyWt,
      postsubset$b_logalpha_Intercept[i]+postsubset$b_logalpha_Log2BodyWt[i]*predout$Log2BodyWt
    )
    lines(predout$Log2BodyWt, samp_lifespan, col = "#3D6FB605", lwd=10);
  }
  dev.off();
}


#posterior predictive check
pp_check(fit_gomp_phylo, resp = "beta", ndraws = 100)
pp_check(fit_gomp_phylo, resp = "logalpha", ndraws = 100)

#forest plot of fit_gomp
library(bayesplot)

# Extract posterior samples for key parameters
posterior_samples <- as.data.frame(fit_gomp_phylo)


# Forest plot of all main effects
mcmc_plot(fit_gomp_phylo, type="hist",
          variable = c("b_beta_Intercept", "b_beta_Log2BodyWt", 
                       "b_logalpha_Intercept", "b_logalpha_Log2BodyWt"),
          prob = 0.9) +
  ggtitle("Histogram of Body Weight Effects")

# Summary table of the parameters
summary(fit_gomp_phylo)

