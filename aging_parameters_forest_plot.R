colors <- c(
            "LifeTable_brms" = "blue", 
            "LifeTable_brms_phylo" = "purple", 
            "LifeTable_lm" = "black", 
            "BrownianVCV"="orange",
            "DiagVCV"="red"
           )



ALLFITS = data.frame()

if ( "LMFIT" %in% ls()) { ALLFITS=bind_rows(ALLFITS,LMFIT) }
  
if ( "brms_alpha" %in% ls()) {
  ALLFITS=bind_rows(ALLFITS,brms_alpha_phylo,brms_beta_phylo)
}

if("STAN_DIAG" %in% ls()) {
   print("Including STAN estimates")
   stan_alpha = data.frame(species=STAN_DIAG$species,  term="alpha0",  estimate=STAN_DIAG$logA, std.error=STAN_DIAG$alphaStdDev, 
                           model="DiagVCV");
   stan_beta = data.frame(species=STAN_DIAG$species,  term="beta0",  estimate=STAN_DIAG$beta, std.error=STAN_DIAG$betaStdDev, 
                          model="DiagVCV")

    ALLFITS = bind_rows(ALLFITS,stan_alpha,stan_beta);
    pdf_filename="primate_forest_plots_stan.pdf";
}
 
if("STAN_BROWN" %in% ls()) {
   print("Including STAN estimates")
   stan_alpha = data.frame(species=STAN_BROWN$species, 
                           term="alpha0",  estimate=STAN_BROWN$logA, std.error=STAN_BROWN$alphaStdDev, 
                           model="BrownianVCV");
   stan_beta = data.frame(species=STAN_BROWN$species,  term="beta0",  estimate=STAN_BROWN$beta, std.error=STAN_BROWN$betaStdDev, 
                          model="BrownianVCV")

    ALLFITS = bind_rows(ALLFITS,stan_alpha,stan_beta);
    pdf_filename="primate_forest_plots_stan.pdf";
} 
## normalize names
ALLFITS$term[ ALLFITS$term == "(Intercept)" ] <- "alpha0"
ALLFITS$term[ ALLFITS$term == "time" ] <- "beta0"
ALLFITS$term[ ALLFITS$term == "INDICESMale" ] <- "alpha_Male"

TREE_ORDER=data.frame(species=paste(primate_tree$tip.label));
TREE_ORDER$tree_order =1:nrow(TREE_ORDER);
tmp=ALLFITS %>% left_join(TREE_ORDER) %>% arrange(desc(tree_order));

p0=ggforestplot::forestplot(
  df = tmp %>% filter(term=="beta0") %>% as.data.frame(),
  name=species,
  estimate = estimate,
  colour = model,
  se = std.error,
  psignif = 0.00001,
  title= "Aging Rate (beta)",
);
p0=p0+scale_color_manual(values = colors) + theme_bw(base_size = 14) +
  geom_vline(xintercept = 0.12 )


p1=ggforestplot::forestplot(
  df = tmp %>% filter(term == "alpha0") %>% as.data.frame(),
  name=species,
  estimate = estimate,
  colour = model,
  se = std.error,
  psignif = 0.00001,
  title = "Baseline Fitness (alpha)",
);
p1=p1+scale_color_manual(values = colors) + theme_bw(base_size = 14) +
  geom_vline(xintercept = -4 ) + 
  theme(legend.position="none"); #turn off legend for one of the plots


pdf_filename = "plots/primate_forest_plots.pdf";
pdf(pdf_filename,width=16,height=8)
p3=ggarrange(p1,p0)
print(p3)
dev.off();

library(ggrepel)
tmp=ALLFITS %>% filter(term %in% c("alpha0","beta0")) %>% 
  reshape2::dcast(species+model~term, value.var = "estimate", fun.aggregate = mean)

t1=reshape2::dcast(ALLFITS %>% filter(term == "beta0"), species ~ term+model, value.var = "estimate", fun.aggregate = median) 
t2=reshape2::dcast(ALLFITS %>% filter(term == "beta0"), species ~ term+model, value.var = "std.error",fun.aggregate = median)
t3=reshape2::dcast(ALLFITS %>% filter(term == "alpha0"), species ~ term+model, value.var = "estimate",fun.aggregate = median)
t4=reshape2::dcast(ALLFITS %>% filter(term == "alpha0"), species ~ term+model, value.var = "std.error",fun.aggregate = median)

##paris plot
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
reg <- function(x, y, col) abline(lm(y~x), col=col) 
panel.lm =  function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
  cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...)  {
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) reg(x[ok], y[ok], col.smooth)
}

 pdf("plots/method_correlations_aging_rates.pdf",width=8,height=8)
 par(mfrow=c(1,2),cex.main=1.3,cex.lab=1.1,pch=19)
 colnames(t1)=c("species","Bayesian"," Bayesian VCV","Lifetable VCV","Lifetable");
 pairs(t1[,c(2:5)],lower.panel = panel.lm, upper.panel = panel.cor,lwd=2, pch=19, 
       main="A. Comparison of Aging Rate Estimates")
 dev.off()
 pdf("plots/method_correlations_baseline_hazards.pdf",width=8,height=8)
 colnames(t3)=c("species","Bayesian"," Bayesian VCV","Lifetable VCV","Lifetable");
 pairs(t3[,c(2:5)],lower.panel = panel.lm, upper.panel = panel.cor,lwd=2, pch=19, 
       main="B. Comparison of Baseline Hazard Estimates")
 dev.off()

 