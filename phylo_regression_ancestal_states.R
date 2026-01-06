# Load necessary libraries
library(ggplot2)
library(dplyr)
library(phytools)
library(nlme)
library(phylolm)
library(flexsurv)
library(latex2exp)
library(RColorBrewer)
library(xtable)

#Bootstrap plots
set.seed(20250101)

# Create a data frame to store the order of species in the tree
TREE_ORDER <- data.frame(species = paste(primate_tree$tip.label))
TREE_ORDER$tree_order <- 1:nrow(TREE_ORDER)

# Join the tree order information to the SPECIES_SUMMARY data frame
tmp <- SPECIES_SUMMARY %>% left_join(TREE_ORDER, by = "species")
SPECIES_SUMMARY$tree_order <- tmp$tree_order

# Create a PDF of the sex maturity distribution
pdf("plots/sex_maturity_distribution.pdf", width = 7.5, height = 5.5)
ggplot(SPECIES_SUMMARY, aes(x = 0, reorder(species, tree_order))) +
  geom_col(aes(x = maxLifespan, fill = "Max"), color = "black", width = 0.5) +
  geom_col(aes(x = quantile.50, fill = "Median"), color = "black", width = 0.5) +
  geom_col(aes(x = SexMature, fill = "Sex Maturity"), color = "black", width = 0.5) +
  scale_fill_manual(values = c("Max" = "#E81313", "Median" = "#FCCE14", "Sex Maturity" = "#1DE213")) +
  xlab("Age") +
  ylab("") +
  theme_bw(base_size = 14) +
  theme(legend.position = c(0.85, 0.85), legend.title = element_blank())
dev.off()

# Summarize alpha0 and beta0 estimates
alpha0 <- ALLFITS %>% 
  filter( model %in% c("BrownianVCV", "LifeTable_brms_phylo") ) %>% 
  filter(term == "alpha0") %>%
  group_by(species) %>%
  dplyr::summarise(alpha = mean(estimate, na.rm = TRUE), alpha_sd=mean(std.error, na.rm=TRUE))


beta0 <- ALLFITS %>%
  filter( model %in% c("BrownianVCV", "LifeTable_brms_phylo") ) %>% 
  filter(term == "beta0") %>%
  group_by(species) %>%
  dplyr::summarise(beta = mean(estimate, na.rm = TRUE), beta_sd=mean(std.error, na.rm=TRUE))

# Combine species summary data with alpha0 and beta0 estimates
SPECIES_SUMMARY_AB <- SPECIES_SUMMARY %>%
  left_join(alpha0, by = "species") %>%
  left_join(beta0, by = "species") %>%
  dplyr::select(
    species, maxLifespan, quantile.50, medianBodyWt, medianBodyWtFemale, medianBodyWtMale, alpha, alpha_sd, beta, beta_sd
  )

# Correct body weight data for specific species
SPECIES_SUMMARY_AB$medianBodyWt[SPECIES_SUMMARY_AB$species == "Owl Monkey"] <- 1.2
SPECIES_SUMMARY_AB$medianBodyWtFemale[SPECIES_SUMMARY_AB$species == "Owl Monkey"] <- 1.2
SPECIES_SUMMARY_AB$medianBodyWtMale[SPECIES_SUMMARY_AB$species == "Owl Monkey"] <- 1.2

SPECIES_SUMMARY_AB$medianBodyWt[SPECIES_SUMMARY_AB$species == "Black howler"] <- 4.4
SPECIES_SUMMARY_AB$medianBodyWtFemale[SPECIES_SUMMARY_AB$species == "Black howler"] <- 5.5
SPECIES_SUMMARY_AB$medianBodyWtMale[SPECIES_SUMMARY_AB$species == "Black howler"] <- 6.7

SPECIES_SUMMARY_AB$medianBodyWt[SPECIES_SUMMARY_AB$species == "Cottontop Tamarin"] <- 0.404
SPECIES_SUMMARY_AB$medianBodyWtFemale[SPECIES_SUMMARY_AB$species == "Cottontop Tamarin"] <- 0.411
SPECIES_SUMMARY_AB$medianBodyWtMale[SPECIES_SUMMARY_AB$species == "Cottontop Tamarin"] <- 0.418

SPECIES_SUMMARY_AB$medianBodyWtMale[SPECIES_SUMMARY_AB$species == "Squirrel Monkey"] <- 0.79
SPECIES_SUMMARY_AB$medianBodyWtFemale[SPECIES_SUMMARY_AB$species == "Squirrel Monkey"] <- 0.68
SPECIES_SUMMARY_AB$medianBodyWtMale[SPECIES_SUMMARY_AB$species == "Squirrel Monkey"] <- 0.899

# Convert columns to numeric type
SPECIES_SUMMARY_AB$beta <- as.double(SPECIES_SUMMARY_AB$beta)
SPECIES_SUMMARY_AB$alpha <- as.double(SPECIES_SUMMARY_AB$alpha)
SPECIES_SUMMARY_AB$medianBodyWt <- as.double(SPECIES_SUMMARY_AB$medianBodyWt)

# Remove rows with missing quantile.50 values
SPECIES_SUMMARY_AB <- SPECIES_SUMMARY_AB %>% filter(!is.na(quantile.50))

# Merge summary statistics with common scientific names
X <- common_sci_names %>%
  dplyr::left_join(SPECIES_SUMMARY_AB, by = c("common_name" = "species"))

# Set row names and convert columns to numeric
rownames(X) <- X$common_name
X$quantile.50 <- as.numeric(X$quantile.50)
X$maxLifespan <- as.numeric(X$maxLifespan)

# Calculate log2 body weight
X$Log2BodyWt <- log2(X$medianBodyWt)
X$Log2BodyWtFemale <- log2(X$medianBodyWtFemale)
X$Log2BodyWtMale <- log2(X$medianBodyWtMale)

# Select relevant columns and remove rows with missing quantile.50 values
X <- X[, c("sci_name", "maxLifespan", "quantile.50", "Log2BodyWt", "Log2BodyWtFemale", "Log2BodyWtMale", "alpha", "alpha_sd", "beta", "beta_sd", "SexMature")]
X <- X[!is.na(X$quantile.50), ]

# Subset the primate tree
primate_tree_subset <- keep.tip(primate_tree, tip = rownames(X))
primate_tree_subset$node.label <- sprintf("%dn", seq(1, length(primate_tree_subset$node.label)) + 10)


# Impute missing data and calculate MRDR
XIMP <- X
XIMP$mrdr <-   log(2) / XIMP$beta
XIMP$mrdr_lb = log(2) / (XIMP$beta + XIMP$beta_sd)
XIMP$mrdr_ub=  log(2) / (XIMP$beta - XIMP$beta_sd)
XIMP$mrdr_sd=  apply(XIMP[,c("mrdr_ub","mrdr_lb")],FUN=sd, MAR=1)


XIMP$species <- rownames(XIMP)
XIMP$logMaxLifespan <- log2(XIMP$maxLifespan)
XIMP$logQuantile50 <- log2(XIMP$quantile.50)

# Summary table for export
XIMP[, c("sci_name", "SexMature", "Log2BodyWt", "alpha", "alpha_sd", "beta", "beta_sd", "mrdr", "mrdr_sd", "quantile.50", "maxLifespan")] %>% clipr::write_clip()
XIMP[, c("sci_name", "SexMature", "Log2BodyWt", "alpha", "alpha_sd", "beta", "beta_sd", "mrdr", "mrdr_sd", "quantile.50", "maxLifespan")] %>% View()
xtable::xtable(XIMP[, c("sci_name", "SexMature", "Log2BodyWt", "alpha", "alpha_sd", "beta", "beta_sd", "mrdr", "mrdr_sd", "quantile.50", "maxLifespan")])

# Phylogenetic signal tests
phylosig(primate_tree_subset, XIMP$Log2BodyWt, method = "K", test = TRUE, nsim = 10000)
phylosig(primate_tree_subset, XIMP$maxLifespan, method = "K", test = TRUE, nsim = 10000)
phylosig(primate_tree_subset, XIMP$quantile.50, method = "K", test = TRUE, nsim = 10000)
phylosig(primate_tree_subset, XIMP$alpha, method = "K", test = TRUE, nsim = 10000)
phylosig(primate_tree_subset, XIMP$mrdr, method = "K", test = TRUE, nsim = 10000)
phylosig(primate_tree_subset, XIMP$beta, method = "K", test = TRUE, nsim = 10000)
phylosig(primate_tree_subset, XIMP$SexMature, method = "K", test = TRUE, nsim = 10000)

phylosig(primate_tree_subset, XIMP$Log2BodyWt, method = "lambda", test = TRUE)
phylosig(primate_tree_subset, XIMP$maxLifespan, method = "lambda", test = TRUE)
phylosig(primate_tree_subset, XIMP$alpha, method = "lambda", test = TRUE)
phylosig(primate_tree_subset, XIMP$quantile.50, method = "lambda", test = TRUE)
phylosig(primate_tree_subset, XIMP$mrdr, method = "lambda", test = TRUE)
phylosig(primate_tree_subset, XIMP$beta, method = "lambda", test = TRUE)
phylosig(primate_tree_subset, XIMP$SexMature, method = "lambda", test = TRUE)

# Ancestral State Reconstruction (MCMC)
par(mfrow = c(1, 1))
COLORS <- c("#E81313", "#FCCE14", "gray", "#6DB1FF", "#1071E5")
Lab.palette <- colorRampPalette(RColorBrewer::brewer.pal(10, "Spectral"))
Lab.palette_inverse <- colorRampPalette(rev(RColorBrewer::brewer.pal(10, "Spectral")))

pdf("plots/primate_tree_out.pdf", width = 7.5, height = 9.5)
par(mfrow = c(1, 1))
for (featureName in c("maxLifespan", "quantile.50", "Log2BodyWt", "alpha", "beta", "mrdr")) {
  roundDigits <- 2
  selectedPallete <- Lab.palette
  if (featureName == "quantile.50") { roundDigits <- 0 }
  if (featureName == "maxLifespan") { roundDigits <- 0 }
  if (featureName == "beta") { selectedPallete <- Lab.palette_inverse }
  if (featureName == "alpha") { selectedPallete <- Lab.palette_inverse }

  Xfeature <- XIMP[, featureName]
  names(Xfeature) <- rownames(XIMP)
  limits <- as.numeric(quantile(Xfeature, c(0.05, 0.95)))
  
  
  #Max Likelihood Ancestral State Reconstruction
  
  corfunc <- corBrownian(1, phy = primate_tree_subset)
  acefit_BM <- ace(Xfeature, primate_tree_subset, corStruct = corfunc, method="ML")
  
  corfunc <- corPagel(0.5, phy = primate_tree_subset)
  acefit_Pagel <- ace(Xfeature, primate_tree_subset, corStruct = corfunc, method="ML")
 
  
  
  #Bayesian MCMC Ancestral State Reconstruction
  mcmc.trait <- anc.Bayes(primate_tree_subset, Xfeature, ngen = 200000)
  mcmc.trait$mean <- apply( mcmc.trait$mcmc, 2, mean)[3:(primate_tree_subset$Nnode + 2)]
  
  data.frame(acefit_BM$ace, acefit_BM$CI95, acefit_Pagel$ace, acefit_Pagel$CI95,mcmc.trait$mean) %>% clipr::write_clip()
  

  cmapout <- contMap(
    primate_tree_subset,
    Xfeature,
    plot = FALSE,
    lim = limits, lwd = 6,
    node.numbers = TRUE,
    outline = TRUE,
    leg.txt = featureName,
    fsize = 0.9
  )

  plot(setMap(cmapout, colors = selectedPallete(100)), lwd = 6, leg.txt = featureName)
  tiplabels(round(Xfeature, roundDigits), adj = c(+0.8, +0.5), bg = "white", cex = 0.5)
  #nodelabels(round(acefit$ace, roundDigits), thermo = acefit$lik.anc, cex = 0.5, bg = "white")
  nodelabels(round(mcmc.trait$mean, roundDigits), cex = 0.5, bg = "white")
}
dev.off()


# MCMC Ancestral State Reconstruction
if (TRUE) {
  trait <- XIMP[, "SexMature"]
  names(trait) <- rownames(XIMP)
  mcmc.sexmature <- anc.Bayes(primate_tree_subset, trait, ngen = 200000)
  mcmcfit.mean.sexmature <- apply(mcmc.sexmature$mcmc, 2, mean)[3:(primate_tree_subset$Nnode + 2)]
  mcmcfit.mean.sexmature_sd <- apply(mcmc.sexmature$mcmc, 2, sd)[3:(primate_tree_subset$Nnode + 2)]

  trait <- XIMP[, "quantile.50"]
  names(trait) <- rownames(XIMP)
  mcmc.q50 <- anc.Bayes(primate_tree_subset, trait, ngen = 200000)
  mcmcfit.mean.q50 <- apply(mcmc.q50$mcmc, 2, mean)[3:(primate_tree_subset$Nnode + 2)]
  mcmcfit.mean.q50_sd <- apply(mcmc.q50$mcmc, 2, sd)[3:(primate_tree_subset$Nnode + 2)]

  trait <- XIMP[, "mrdr"]
  names(trait) <- rownames(XIMP)
  mcmc.mrdr <- anc.Bayes(primate_tree_subset, trait, ngen = 200000)
  mcmcfit.mean.mrdr <- apply(mcmc.mrdr$mcmc, 2, mean)[3:(primate_tree_subset$Nnode + 2)]
  mcmcfit.mean.mrdr_sd <- apply(mcmc.mrdr$mcmc, 2, sd)[3:(primate_tree_subset$Nnode + 2)]
  

  trait <- XIMP[, "alpha"]
  names(trait) <- rownames(XIMP)
  mcmc.alpha <- anc.Bayes(primate_tree_subset, trait, ngen = 200000)
  mcmcfit.mean.alpha <- apply(mcmc.alpha$mcmc, 2, mean)[3:(primate_tree_subset$Nnode + 2)]
  mcmcfit.mean.alpha_sd <- apply(mcmc.alpha$mcmc, 2, sd)[3:(primate_tree_subset$Nnode + 2)]

  trait <- XIMP[, "beta"]
  names(trait) <- rownames(XIMP)
  mcmc.beta <- anc.Bayes(primate_tree_subset, trait, ngen = 200000)
  mcmcfit.mean.beta <- apply(mcmc.beta$mcmc, 2, mean)[3:(primate_tree_subset$Nnode + 2)]
  mcmcfit.mean.beta_sd <- apply(mcmc.beta$mcmc, 2, sd)[3:(primate_tree_subset$Nnode + 2)]

  trait <- XIMP[, "maxLifespan"]
  names(trait) <- rownames(XIMP)
  mcmc.maxLifespan <- anc.Bayes(primate_tree_subset, trait, ngen = 200000)
  mcmcfit.mean.maxLifespan <- apply(mcmc.maxLifespan$mcmc, 2, mean)[3:(primate_tree_subset$Nnode + 2)]
  mcmcfit.mean.maxLifespan_sd <- apply(mcmc.maxLifespan$mcmc, 2, sd)[3:(primate_tree_subset$Nnode + 2)]

  projQ50 <- qgompertz(0.5, shape = mcmcfit.mean.beta, rate = exp(mcmcfit.mean.alpha))

  data.frame(
    node = primate_tree_subset$node.label,
    alpha = mcmcfit.mean.alpha,
    alpha_sd = mcmcfit.mean.alpha_sd,
    beta = mcmcfit.mean.beta,
    beta_stddev = mcmcfit.mean.beta_sd,
    mrdr = mcmcfit.mean.mrdr,
    mrdr_stddev = mcmcfit.mean.mrdr_sd,
    medianLifespan = mcmcfit.mean.q50,
    medianLifespan_sd = mcmcfit.mean.q50_sd,
    maxLifespan = mcmcfit.mean.maxLifespan,
    maxLifAespan_stddev = mcmcfit.mean.maxLifespan_sd
  ) %>% write.table("outputs/ancestral_states.tsv", sep = "\t")
}


if (1) {
  pdf("plots/ancestral_states_recon.pdf", width = 12, height = 7.5)
  par(mfrow = c(1, 2))

  trait <- XIMP[, "mrdr"]
  names(trait) <- rownames(XIMP)
  limits <- as.numeric(quantile(trait, c(0.05, 0.95)))
  cmapout <- contMap(primate_tree_subset, trait, lims = limits, plot = FALSE, CE=TRUE)
  nodelabel <- paste0(round(mcmcfit.mean.mrdr, 1))
  rounding <- 1
  
  plot(setMap(cmapout, colors = Lab.palette(100)), lwd = 6, sig = 1, node.numbers = TRUE, outline = TRUE, font = 1, offset = 1.5, cex = 0.8, ftype = "reg",
       fsize = 0.9, edge.width = 8, leg.txt = "MRDT (yrs.)")

  nodelabels(nodelabel, cex = 0.7, bg = "#FFFFFF", frame = "rect")
  tiplabels(round(trait, rounding), adj = c(-0.4, +0.5), bg = "#FFFFFF00", cex = 0.7, frame = "rect")

  trait <- XIMP[, "alpha"]
  names(trait) <- rownames(XIMP)
  limits <- as.numeric(quantile(trait, c(0.05, 0.95)))
  cmapout <- contMap(primate_tree_subset, trait, lims = limits, plot = FALSE, CE=TRUE)

  plot(setMap(cmapout, colors = Lab.palette_inverse(100)), lwd = 6, sig = 1, node.numbers = TRUE, outline = TRUE, font = 1, offset = 2, cex = 0.8, ftype = "reg", fsize = 0.9,
       edge.width = 8, leg.txt = "Baseline Hazard")

  nodelabel <- paste0(round(mcmcfit.mean.alpha, 1))
  nodelabels(nodelabel, cex = 0.7, bg = "#FFFFFF", frame = "rect")
  tiplabels(round(trait, 1), adj = c(-0.4, +0.5), bg = "#FFFFFF00", cex = 0.7, frame = "rect")
  dev.off()
}


if (1) {
  pdf("plots/ancestral_states_lifespan.pdf", width = 12, height = 8)
  par(mfrow = c(1, 2))
  rounding <- 0

  trait <- XIMP[, "quantile.50"]
  names(trait) <- rownames(XIMP)
  limits <- as.numeric(quantile(trait, c(0.05, 0.95)))
  cmapout <- contMap(primate_tree_subset, trait, lims = limits, plot = FALSE, CE=TRUE)
  
  plot(setMap(cmapout, colors = Lab.palette(100)), lwd = 6, sig = 1, node.numbers = TRUE, outline = TRUE, font = 1, offset = 2, cex = 0.8, ftype = "reg", fsize = 0.9,
       edge.width = 8, leg.txt = "Median Lifespan (yrs.)")
  
  nodelabel <- paste0(round(mcmcfit.mean.q50, 0))
  nodelabels(nodelabel, cex = 0.7, bg = "#FFFFFF", frame = "rect")
  tiplabels(round(trait, 1), adj = c(-0.4, +0.5), bg = "#FFFFFF00", cex = 0.7, frame = "rect")
  
  
  trait <- XIMP[, "maxLifespan"]
  names(trait) <- rownames(XIMP)
  limits <- as.numeric(quantile(trait, c(0.05, 0.95)))
  cmapout <- contMap(primate_tree_subset, trait, lims = limits, plot = FALSE, CE=TRUE)
  nodelabel <- paste0(round(mcmcfit.mean.maxLifespan, 0))

  
  plot(setMap(cmapout, colors = Lab.palette(100)), lwd = 6, sig = 1, node.numbers = TRUE, outline = TRUE, font = 1, offset = 1.5, cex = 0.8, ftype = "reg",
       fsize = 0.9, edge.width = 8, leg.txt = "Maximum Lifespan (yrs.)")
  
  nodelabels(nodelabel, cex = 0.7, bg = "#FFFFFF", frame = "rect")
  tiplabels(round(trait, rounding), adj = c(-0.4, +0.5), bg = "#FFFFFF00", cex = 0.7, frame = "rect")
  
  dev.off()
}



# Phylogenetic Regressions

bootStrapLM <- function(formula_str, df, nboot = 100, ncoef = 2) {
  COEF <- matrix(nrow = nboot, ncol = ncoef)
  for (i in 1:nboot) {
    dfsubset <- df %>% dplyr::sample_frac(0.7)
    lmfit <- lm(as.formula(formula_str), dfsubset)
    COEF[i, ] <- coef(lmfit)
  }
  return(COEF)
}

bootStrapGLS <- function(formula_str, df, nboot = 100, ncoef = 2) {
  COEF <- matrix(nrow = nboot, ncol = ncoef)
  for (i in 1:nboot) {
    dfsubset <- df %>% dplyr::sample_frac(0.8)
    tmptree <- keep.tip(primate_tree_subset, tip = rownames(dfsubset))
    dfsubset$species <- rownames(dfsubset)
    corMatrix <- corBrownian(1, phy = tmptree, form = ~species)
    lmfit <- gls(as.formula(formula_str), dfsubset, correlation = corMatrix)
    COEF[i, ] <- coef(lmfit)
  }
  return(COEF)
}

bootStrapPhylom <- function(formula_str, df, nboot = 100, ncoef = 2) {
  COEF <- matrix(nrow = nboot, ncol = ncoef)
  for (i in 1:nboot) {
    dfsubset <- df %>% dplyr::sample_frac(0.8)
    tmptree <- keep.tip(primate_tree_subset, tip = rownames(dfsubset))
    dfsubset$species <- rownames(dfsubset)
    pfit <- phylolm::phylolm(formula_str, data = dfsubset, phy = tmptree, model = "BM")
    COEF[i, ] <- coef(pfit)
  }
  return(COEF)
}

# Model Prediction Evaluation

if (TRUE) {
  calcR2 <- function(yres, y, dfr = 37, dft = 39) {
    yres <- sum(yres^2)
    ytol <- sum((y - mean(y))^2)
    1 - ((yres / dfr) / (ytol / dft))
  }

  XIMP$maxAdultLifespan <- XIMP$maxLifespan - XIMP$SexMature
  XIMP$medAdultLifespan <- XIMP$quantile.50 - XIMP$SexMature

  primate_tree_subset$root.edge <- 0
  pdf("plots/model_evaluation.pdf", width = 10, height = 5)

  BOOT <- 100
  par(mfrow = c(1, 2), cex.main = 1.1, cex = 1.0, cex.lab = 1.05)

  XIMP$predQ99 <- qgompertz(0.99, shape = XIMP$beta, rate = exp(XIMP$alpha))
  plot(predQ99 ~ maxAdultLifespan, data = XIMP, pch = 19,
       main = "A. Maximum Lifespan",
       xlab = "Maximum Adult Lifespan",
       ylab = "Predicted Max. Adult Lifespan")
  abline(c(0, 1), lty = 2, lwd = 3)
  r2 <- calcR2(XIMP$predQ99 - XIMP$maxAdultLifespan, XIMP$maxAdultLifespan)
  legend("bottomright", latex2exp::TeX(sprintf("$R^2$=%3.2f", round(r2, 2))))

  XIMP$predQ50 <- qgompertz(0.5, shape = XIMP$beta, rate = exp(XIMP$alpha))
  plot(predQ50 ~ medAdultLifespan, data = XIMP, pch = 19,
       main = "B. Median Lifespan",
       xlab = "Median Adult Lifespan",
       ylab = "Predicted Median Adult Lifespan")
  abline(c(0, 1), lty = 2, lwd = 3)
  r2 <- calcR2(XIMP$predQ50 - XIMP$medAdultLifespan, XIMP$medAdultLifespan)
  legend("bottomright", latex2exp::TeX(sprintf("$R^2$=%3.2f", round(r2, 2))))

  dev.off()
}

if (1) {
  pdf("plots/fixed_alpha_beta.pdf", width = 10, height = 10)
  par(mfrow = c(2, 2), cex.main = 1.1, cex = 1.0, cex.lab = 1.05)

  XIMP$predQ50 <- qgompertz(0.5, shape = 0.13, rate = exp(XIMP$alpha))
  calcR2(XIMP$predQ50, XIMP$maxAdultLifespan)
  plot(predQ50 ~ medAdultLifespan, data = XIMP, pch = 19,
       main = expression(paste("A. Median Lifespan: Fixed ", beta, "=0.13")),
       xlab = "Median Adult Lifespan",
       ylab = "Predicted Median Adult Lifespan")
  abline(c(0, 1), lty = 2, lwd = 3)
  r2 <- calcR2(XIMP$predQ50 - XIMP$medAdultLifespan, XIMP$medAdultLifespan)
  legend("bottomright", latex2exp::TeX(sprintf("$R^2$=%3.2f", round(r2, 2))))


  XIMP$predQ99 <- qgompertz(0.99, shape = 0.13, rate = exp(XIMP$alpha))
  calcR2(XIMP$predQ99, XIMP$maxAdultLifespan)
  plot(predQ99 ~ maxAdultLifespan, data = XIMP, pch = 19,
       main = expression(paste("B. Max Lifespan: Fixed ", beta, "=0.13")),
       xlab = "Maximum Adult Lifespan",
       ylab = "Predicted Maximum Adult Lifespan")
  abline(c(0, 1), lty = 2, lwd = 3)
  r2 <- calcR2(XIMP$predQ99 - XIMP$maxAdultLifespan, XIMP$maxAdultLifespan)
  legend("bottomright", latex2exp::TeX(sprintf("$R^2$=%3.2f", round(r2, 2))))



  XIMP$predQ50 <- qgompertz(0.5, shape = XIMP$beta, rate = exp(-4.3))
  plot(predQ50 ~ medAdultLifespan, data = XIMP, pch = 19,
       main = expression(paste("C. Median Lifespan: Fixed ", ln(alpha), "=-4.3")),
       xlab = "Median Adult Lifespan",
       ylab = "Predicted Median Adult Lifespan")
  abline(c(0, 1), lty = 2, lwd = 3)
  r2 <- calcR2(XIMP$predQ50 - XIMP$medAdultLifespan, XIMP$medAdultLifespan)
  legend("bottomright", latex2exp::TeX(sprintf("$R^2$=%3.2f", round(r2, 2))))

  XIMP$predQ99 <- qgompertz(0.99, shape = XIMP$beta, rate = exp(-4.3))
  plot(predQ99 ~ maxAdultLifespan, data = XIMP, pch = 19,
       main = expression(paste("D. Max Lifespan: Fixed ", ln(alpha), "=-4.3")),
       xlab = "Maximum Adult Lifespan",
       ylab = "Predicted Maximum Adult Lifespan")
  abline(c(0, 1), lty = 2, lwd = 3)
  r2 <- calcR2(XIMP$predQ99 - XIMP$maxAdultLifespan, XIMP$maxAdultLifespan)
  legend("bottomright", latex2exp::TeX(sprintf("$R^2$=%3.2f", round(r2, 2))))

  dev.off()
}

# Alpha Beta Correlation
pfitBeta <- phylolm::phylolm(beta ~ alpha,
  data = XIMP,
  phy = primate_tree_subset,
  model = "BM",
  boot = 1000
)
summary(pfitBeta)


doPhyloRegression <- function(yvar, xvar, ylab, xlab, title, modeltype="BM", showUndjusted=FALSE, legendpos="bottomright") {
  formula_str <- as.formula(paste0(yvar, "~", xvar))
  plot(formula_str, data = XIMP, pch = 19, xlab = xlab, ylab = ylab, main = title)

  pfit <- phylolm::phylolm(formula_str, data = XIMP, phy = primate_tree_subset, model = modeltype, boot = 1000)
  pfit$bootstrap <- bootStrapPhylom(formula_str, XIMP, nboot = BOOT)
  bootmean=colMeans(pfit$bootstrap)
  for (i in 1:BOOT) {
    abline(coef = pfit$bootstrap[i, 1:2], col = "#3D6FB605", lwd = 10)
  }
  points(formula_str, data = XIMP, col = "black", pch = 19)
  abline(coef=bootmean, lwd = 5, col = "#3D6FB6")
  print(summary(pfit))
  r2 <- pfit$r.squared
  
  if(showUndjusted) { 
    lmfit <- lm(formula_str, data = XIMP)
    lmfit$bootstrap <- bootStrapLM(formula_str, XIMP, nboot = BOOT)
    for (i in 1:BOOT) {
      abline(coef = lmfit$bootstrap[i, 1:2], col = "#ff1f3d05", lwd = 10)
    }
    abline(coef = colMeans(lmfit$bootstrap), lwd = 3, col = "#ff1f3d", lty = 2)
    tmp <- summary(lmfit)
    lmr2 <- tmp$r.squared
    #legend("bottomright", latex2exp::TeX(sprintf("$R^2$=%3.2f $\\phyR^2$=%3.2f", round(lmr2, 2), round(r2, 2))), border = 0)
  }

  legend(legendpos, latex2exp::TeX(sprintf("$phyR^2$=%3.2f", round(r2, 2))), border = 0)
  addTextLabels(XIMP[, xvar], XIMP[, yvar], rownames(XIMP), cex.label = 0.8, lty=2, col.label = "black")
}

# BodyWeight Correlation
if (1) {
  pdf("plots/body_weight_correlation.pdf", width = 10, height = 10)
  par(mfrow = c(2, 2), cex.lab = 1.2, cex.axis = 1.3, cex.main = 1.4)
  doPhyloRegression("logQuantile50", "Log2BodyWt", ylab = "Log2(Median Lifespan) yrs.", xlab = "Log2(BodyWeight) kg.", title = "A. Median Lifespan vs Body Weight", modeltype="lambda")
  doPhyloRegression("logMaxLifespan", "Log2BodyWt", ylab = "Log2(Maximum Lifespan) yrs.", xlab = "Log2(BodyWeight) kg.", title = "B. Maximum Lifespan vs Body Weight",modeltype="lambda")
  doPhyloRegression("alpha", "Log2BodyWt", ylab = "Baseline Hazard ln(alpha)", xlab = "Log2(BodyWeight) kg.", title = "C. Baseline Hazard vs Body Weight",modeltype="lambda")
  doPhyloRegression("beta", "Log2BodyWt", ylab = "Aging Rate", xlab = "Log2(BodyWeight) kg.", title = "D. Aging Rate vs Body Weight",modeltype="lambda")
  dev.off()
}

## 
if(1) {
pdf("plots/alpha_beta_correlation_models.pdf", width = 10, height = 10)
  par(mfrow = c(2, 2), cex.lab = 1.2, cex.axis = 1.3, cex.main = 1.4)

label1=latex2exp::TeX("Aging Rate $\\beta$")
label2=latex2exp::TeX("Baseline Hazard $ln(\\alpha)$")
doPhyloRegression("beta", "alpha", ylab=label1, xlab=label2, title = "Baseline Hazard vs Aging Rate (Brownian)",modeltype="BM", legendpos = "topleft")
doPhyloRegression("beta", "alpha", ylab=label1, xlab=label2, title = "Baseline Hazard vs Aging Rate (Pagel Lambda)",modeltype="lambda", legendpos = "topleft")
doPhyloRegression("beta", "alpha", ylab=label1, xlab=label2, title = "Baseline Hazard vs Aging Rate (OUfixedRoot)",modeltype="OUfixedRoot", legendpos = "topleft")
doPhyloRegression("beta", "alpha", ylab=label1, xlab=label2, title = "Baseline Hazard vs Aging Rate (OUrandomRoot)",modeltype="OUrandomRoot", legendpos = "topleft")
dev.off()
}

# CORRELATIONS
cor.brown <- corBrownian(1, phy = primate_tree_subset, form = ~species)
pfitBeta <- phylolm::phylolm(beta ~ alpha, data = XIMP, phy = primate_tree_subset, model = "lambda")
summary(pfitBeta)
doPhyloRegression("logMaxLifespan", "alpha", xlab = "Baseline Hazard ln(alpha)", ylab = "max Lifespan", title = "Baseline Hazard vs maxLifespan",modeltype="lambda")
doPhyloRegression("logMaxLifespan", "beta", xlab = "Aging Rate", ylab = "max Lifespan", title = "Aging Rate vs maxLifespan",modeltype="lambda")

pfitMRDR <- phylolm::phylolm(mrdr ~ Log2BodyWt, data = XIMP, phy = primate_tree_subset, model = "lambda")
glsfit <- nlme::gls(mrdr ~ Log2BodyWt, correlation = cor.brown, data = XIMP)
summary(glsfit)
summary(pfitMRDR)

pfitAlpha <- phylolm::phylolm(alpha ~ Log2BodyWt, data = XIMP, phy = primate_tree_subset, model = "BM")
glsfit <- nlme::gls(alpha ~ Log2BodyWt, correlation = cor.brown, data = XIMP)
summary(glsfit)

XIMP$predQ50 <- flexsurv::qgompertz(0.99, shape = predict(pfitBeta), rate = exp(predict(pfitAlpha)))
plot(XIMP$maxLifespan, XIMP$predQ50, pch = 19)
summary(lm(XIMP$maxLifespan ~ XIMP$predQ50))
