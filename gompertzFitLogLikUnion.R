ESTIMATE=FALSE
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
NumMCMCChains=4

##setup model data frame
primate_tree_subset = keep.tip(primate_tree,tip = unique(ALL$species))
dd = ALL %>% filter(species %in% primate_tree_subset$tip.label) %>% 
  filter(lifespan>=sexMature)

dd$lifespan = dd$lifespan - dd$sexMature #adjust t0 for sexual maturity

tree_oder=data.frame(species=paste(primate_tree_subset$tip.label), tree_order=1:(primate_tree_subset$Nnode+1))
dd=dd %>% left_join(tree_oder)

imalanced=dd %>% group_by(species) %>% 
  summarise(n=n(), nFemale=sum(sex=="Female"),  nMale=sum(sex=="Male")) %>% 
  filter(nFemale/nMale < 0.5 | n < 100)

## can't fit..assumes sex==0 
for(sp in unique(imalanced$species)) {
  print(sp);
  dd$sex[dd$species == sp] <- "Female"
}

# stan fits
stanFitSurvivalCurveLogLik =function(modeldf, Sigma) {
  model_data = list(
    N=nrow(modeldf),
    K=length(unique(modeldf$species)),
    species = modeldf$tree_order,
    event=modeldf$event,
    time=modeldf$lifespan,
    sex=ifelse(modeldf$sex == "Female",0,1),
    VCV=Sigma
  );
  
  #starting parameters for alpha and beta
  initf = function() { list("alpha0"=-3, "beta0"=-1) }
  
  #run MCMC
  model <- stan_model(file="gompertzFitLogLikUnion.stan");
  fit <- rstan::sampling(model,
              data = model_data,  
              chains = NumMCMCChains, 
              iter = 500,
              init = initf
              );
  return(fit)
}



#FIT BAYESIAN MODEL

if (ESTIMATE) {
  
  #estimate, no information sharing
  VCV=matrix(0,nrow(VCV),ncol(VCV));  
  diag(VCV)=1; 
  nspp=nrow(VCV); VCV <- VCV/(det(VCV)^(1/nspp))# rescale vcv
  
  #Fit with RSTAN
  stanfitDiag  = stanFitSurvivalCurveLogLik(dd, VCV)
  
  #with information sharing
  VCV=vcv(primate_tree_subset) #BROWNIAN 
  nspp=nrow(VCV); VCV <- VCV/(det(VCV)^(1/nspp)); #rescale vcv
  
  #Fit with RSTAN
  stanfitBrownian  = stanFitSurvivalCurveLogLik(dd, VCV)
  
  #VCV=vcv(corMartins(2,primate_tree_subset))
  #heatmap(VCV)
  
  saveRDS(stanfitDiag,"outputs/gompertzFitLogLikUnion.diag.fit") 
  saveRDS(stanfitBrownian,"outputs/gompertzFitLogLikUnion.brownian.fit")
} else { 
  #read saved data
  stanfitDiag=readRDS("outputs/gompertzFitLogLikUnion.diag.fit")
  stanfitBrownian=readRDS("outputs/gompertzFitLogLikUnion.brownian.fit")
}

#### PLOTS#####

logit=function(x) { log(x/(1-x)); }
invlogit=function(x) {  exp(x)/(1+exp(x)); }

getStanFit=function(stanfit) {
  Z=summary(stanfit)$summary
  ZZ=extract(stanfit) %>% as.data.frame()
  FINAL=tree_oder;
  ##alpha effect
  alphaSpecies=ZZ %>% dplyr::select(dplyr::matches("alphaSpecies\\.",perl=T));
  ## adjust for common intercept if it exists
  #a0=ZZ %>% select(matches("alpha0",perl=T));
  #for(i in 1:ncol(alphaSpecies)) { alphaSpecies[,i] = alphaSpecies[,i]+a0; }
  FINAL$logA = alphaSpecies %>% summarise_all(mean) %>% as.double()
  FINAL$alphaStdDev = alphaSpecies %>% summarise_all(sd) %>% as.double()
  FINAL$alphaSex = Z[grep("alphaSpeciesSex\\[",rownames(Z)),"mean"]
  FINAL$alphaNeff = Z[grep("alphaSpecies\\[",rownames(Z)),"n_eff"]
  
  ### beta effect
  FINAL$beta=ZZ %>% dplyr::select(matches("betaSpecies\\.",perl=T)) %>% 
    mutate_all(invlogit) %>% summarise_all(median) %>% as.double()
  FINAL$betaStdDev = ZZ %>% dplyr::select(matches("betaSpecies\\.",perl=T)) %>% 
    mutate_all(invlogit) %>% summarise_all(sd) %>% as.double()
  FINAL$betaSex = Z[grep("betaSpeciesSex\\[",rownames(Z)),"mean"]
  FINAL$betaNeff = Z[grep("betaSpecies\\[",rownames(Z)),"n_eff"]
  FINAL$mrdr  = log(2)/FINAL$beta;
  rownames(FINAL) = FINAL$species
  return(FINAL);
}

STAN_DIAG=getStanFit(stanfitDiag);
STAN_BROWN=getStanFit(stanfitBrownian);

#FINAL$mu_beta = mu_beta
#mu_beta =exp(Z[grep("mu_beta",rownames(Z)),"mean-all chains"])
#vcv_scale =exp(Z[grep("vcv_scale",rownames(Z)),"mean-all chains"])

pdf("plots/primate_fits_union_tree.pdf",width=10,height=8)
par(mfrow=c(1,1));

FINAL=STAN_DIAG

par(mfrow=c(1,1))
roundDigits=2
for(featureName in c("logA","mrdr")) { 
  Lab.palette <- colorRampPalette(c("red", "gray", "blue"), space = "rgb")
  Xfeature=FINAL[,featureName]; names(Xfeature)= FINAL$species;
  limits=as.numeric(quantile(Xfeature,c(0.05,0.95)))
  acefit=ace(Xfeature,primate_tree_subset,corStruct = VCV,type = "continuous")
  contMap(primate_tree_subset,Xfeature,plot=TRUE,lim=limits,lwd=6,palette=Lab.palette,node.numbers=TRUE,outline=TRUE)
  tiplabels(round(Xfeature,roundDigits),adj = c(+0.8,+0.5),bg="white",cex=0.5)
  nodelabels(round(acefit$ace,roundDigits),thermo = acefit$lik.anc, cex=0.5,bg="white") 
}
dev.off()

pdf("plots/primate_fits_union.pdf",width=11,height=8)
SUBSET=unique(FINAL$species)
for( cohortName in SUBSET) {
  print(cohortName)
  #subset species
  Z = dd %>% dplyr::filter(species == cohortName) %>% arrange(lifespan)
  Z$sex  = ifelse(Z$sex == "Female", 0, 1)
  if(length(unique(Z$sex)) == 1) { Z$sex=0 };

  #FIT BAYESIAN MODEL
  estA=exp(FINAL[cohortName,"logA"])
  estB=FINAL[cohortName,"beta"];
  par(mfrow=c(1,2))
  PERIOD=1
  COLORS=c("black","red");
  haz = function(x,a,b) { a*exp(b*x); };  
  cumhaz = function(x,a,b) { (a/b)*( exp(b*x)-1 ); }
  sur = function(x,a,b) { exp(-cumhaz(x,a,b)); }
    
  maxLifespan=round(max(Z$lifespan)+1)
  alphaSex=FINAL[cohortName,"alphaSex"];
  betaSex=FINAL[cohortName,"betaSex"];
  N=sum(Z$event)
  if(N<20) { next; }
    sfit=survfit(Surv(Z$lifespan,Z$event) ~ Z$sex)
    plot(sfit, conf.int=FALSE,main=cohortName, ylab="Frac Alive", xlab="Age", mark.time=TRUE, cex=1.5, col=COLORS);
    subTitle=sprintf("a=%.5f b=%.4f", estA, estB)
    points(Z$lifespan,sur(Z$lifespan,estA,estB),col=COLORS,type="l",lwd=3)
    points(Z$lifespan,sur(Z$lifespan,estA*exp(alphaSex),estB*exp(betaSex)),col="red",type="l",lwd=3)
    mtext(subTitle,side=3,cex=0.8)
    
    lt = LT %>% filter(species == cohortName);
    #lt <- dd %>% filter(species == cohortName) %>% group_by(species,sex) %>% do(lifetable(.))
    #lt$time = lt$int.start+(lt$int.end-lt$int.start)/2
    
    if (nrow(lt)>0) {
      timedomain=seq(min(lt$int.start), max(lt$int.end),1)
      plot(lt$time,lt$rate,log="y",pch=19, xlab="Age", ylab="Hazard",col=as.factor(lt$sex))
      lines(timedomain,haz(timedomain,estA, estB),col=COLORS[1], lwd=2);
      mtext(sprintf("Ndeath=%d ln(a)= %3.2f DoubleRate= %3.2f", N, log(estA),log(2)/estB), side=3, cex=0.8);
      
      mfits=ALLFITS %>% reshape2::dcast(model +species ~ term, value.var = "estimate",fun.aggregate = mean) %>% 
        filter(species == cohortName)
      
      for( m in 1:nrow(mfits)) { 
        lines(timedomain, haz(timedomain,exp(mfits$alpha0[m]),mfits$beta0[m]), col=m, lwd=2);
      }
      legend("topleft",c(mfits$model),fill=c(1:nrow(mfits)))
    }
}
dev.off()
