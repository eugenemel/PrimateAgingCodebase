
functions {
  real log_h(real T,  real alpha, real beta);
  real H_t(real T,  real alpha, real beta);
  real invLogit(real x);
  
  //return log hazard at time T
  real log_h(real T,  real alpha, real beta) { return log(alpha) + beta*T;  }
  
  //returns cumulative hazard up to time T
  real H_t(real T,  real alpha, real beta){ return (alpha/beta)*( exp(beta*T)-1 ); }
  
  //inverse logit x -> [0,1] range
  real invLogit(real x){  return( exp(x)/(1+exp(x)) ); }
}

data { 
  int<lower=1> N; //number of subjects
  int<lower=1> K; //number of species
  int species[N];
  real time[N];
  real event[N];
  int  sex[N];
  matrix[K,K] VCV;
}

parameters { 
  real<upper=0>alpha0;
  real<lower=0>sigma_alpha0;
  
  real beta0;
  real<lower=0>sigma_beta0;
  
  vector[K] mu_alpha;
  vector<lower=0>[K] sigma_alpha;
  vector<upper=0>[K] alphaSpecies;
  
  vector[K] mu_beta;
  vector[K] betaSpecies; 
  vector<lower=0,upper=1>[K] sigma_beta;
  
  real<lower=0>sigma_alpha_sex;
  real<lower=0>sigma_beta_sex;  
  vector[K] alphaSpeciesSex;
  vector[K] betaSpeciesSex;
  
  real<lower=0,upper=1> vcv_scale_alpha;
  real<lower=0,upper=1> vcv_scale_beta;
  real<lower=0>sigma_vcv_scale_alpha;
  real<lower=0>sigma_vcv_scale_beta;  
}
 
transformed parameters {
  matrix[K,K] alphaVCV=VCV;
  matrix[K,K] betaVCV=VCV;

  for(r in 1:K) { 
     for(k in 1:K) { 
       if(r != k ) { //off diaganal elements
          alphaVCV[r,k]= VCV[r,k]*vcv_scale_alpha;  //Pagel
          betaVCV[r,k]=  VCV[r,k]*vcv_scale_beta;   //Pagel
          //alphaVCV[r,k]= pow(VCV[r,k],vcv_scale_alpha);
          //betaVCV[r,k]=  pow(VCV[r,k],vcv_scale_beta);
       }
   }
  }
  
  print("a=",alphaSpecies);
  print("b=",betaSpecies);
}

model {
  
  for (i in 1:N) {
        int sid=species[i]; //get species for an individual
        //baseline hazard
        real A = exp(alphaSpecies[sid]+alphaSpeciesSex[sid]*sex[i]);
        
        //aging rate
        real B = invLogit(betaSpecies[sid]+betaSpeciesSex[sid]*sex[i]);
        
        //gompertz lpdf
        target += event[i]*log_h(time[i], A, B) - H_t(time[i], A, B);
  }
  
  //target += lkj_corr_cholesky_lpdf(betaVCV | 1);
  //target += lkj_corr_cholesky_lpdf(alphaVCV | 1);
  target += std_normal_lpdf(mu_alpha);
  target += std_normal_lpdf(mu_beta);
    
  //VCV scalling
  sigma_vcv_scale_alpha ~ normal(0,1);
  sigma_vcv_scale_beta ~ normal(0,1);
  vcv_scale_alpha ~ normal(0,sigma_vcv_scale_alpha);
  vcv_scale_beta~ normal(0, sigma_vcv_scale_beta);

  //global fit
  sigma_alpha0 ~ normal(0,5);
  alpha0 ~ normal(0, sigma_alpha0);
 
  //global slope 
  sigma_beta0 ~ normal(0,2);
  beta0 ~  normal(0, sigma_beta0);
  
  //species
  sigma_alpha ~ normal(0,3);
  sigma_beta ~ normal(0,3);
  mu_alpha ~ normal(alpha0, sigma_alpha);
  mu_beta  ~ normal(beta0, sigma_beta);
  betaSpecies ~  multi_normal(mu_beta,betaVCV);
  alphaSpecies  ~  multi_normal(mu_alpha,alphaVCV);
  
  //sex adjustment
  sigma_alpha_sex ~ normal(0,1);
  sigma_beta_sex ~  normal(0,1);
  alphaSpeciesSex ~ double_exponential(0,sigma_alpha_sex);
  betaSpeciesSex ~  double_exponential(0,sigma_beta_sex);
}

//must have blank line below
