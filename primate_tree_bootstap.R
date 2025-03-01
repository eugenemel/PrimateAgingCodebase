# Set up the bootstrap parameters
NBOOT <- 100

# Initialize empty data frames for bootstrap results
ACEBOOTSTRAP <- data.frame()
KMBOOTSTRAP <- data.frame()

# Filter species with available median lifespan data
keepSpecies <- SPECIES_SUMMARY$species[!is.na(SPECIES_SUMMARY$quantile.50)]

# Subset the phylogenetic tree to keep only relevant species
primate_tree_subset <- keep.tip(primate_tree, tip = keepSpecies)
primate_tree_subset$node.label <- sprintf("%dn", seq_along(primate_tree_subset$node.label) + 10)

# Visualize the tree
plot(primate_tree_subset, show.node.label = TRUE)

# Compute phylogenetic distances using Brownian motion model
cor.brown <- corBrownian(1, phy = primate_tree_subset)

# Run bootstrap procedure
for (i in 1:NBOOT) {
  print(paste("Bootstrap iteration:", i))
  
  # Random subset of lifespan data (80% sampling)
  tmp <- ALL %>% 
    filter(species %in% keepSpecies) %>% 
    group_by(species) %>% 
    sample_frac(0.8)
  
  # Kaplan-Meier survival analysis
  kmfit <- survfit(Surv(lifespan, event) ~ species, data = tmp)
  
  # Extract median and quantile information
  medianKMSurv <- as.data.frame(quantile(kmfit, probs = c(0.5, 0.95))) %>%  
    dplyr::select(quantile.50, quantile.95)
  medianKMSurv$species <- gsub("species=", "", rownames(medianKMSurv))
  
  # Calculate maximum lifespan for each species
  maxLifespan <- tmp %>% 
    group_by(species) %>% 
    summarise(maxLifespan = max(lifespan))
  
  # Join the datasets
  medianKMSurv <- medianKMSurv %>% left_join(maxLifespan)
  
  # Process maximum lifespan ancestral states
  X <- medianKMSurv$maxLifespan
  if (any(is.na(X))) next
  
  names(X) <- medianKMSurv$species
  acefit <- ace(X, primate_tree_subset, corStruct = cor.brown, method = "GLS")
  
  # Save maximum lifespan results
  ACEBOOTSTRAP <- rbind(
    ACEBOOTSTRAP,
    data.frame(node = names(acefit$ace), ace = acefit$ace, iter = i, estimate = "maximum")
  )
  
  KMBOOTSTRAP <- rbind(
    KMBOOTSTRAP,
    data.frame(node = names(X), X = X, iter = i, estimate = "maximum")
  )
  
  # Process median lifespan ancestral states
  X <- medianKMSurv$quantile.50
  if (any(is.na(X))) next
  
  names(X) <- medianKMSurv$species
  acefit <- ace(X, primate_tree_subset, corStruct = cor.brown, method = "GLS")
  
  # Save median lifespan results
  ACEBOOTSTRAP <- rbind(
    ACEBOOTSTRAP,
    data.frame(node = names(acefit$ace), ace = acefit$ace, iter = i, estimate = "median")
  )
  
  KMBOOTSTRAP <- rbind(
    KMBOOTSTRAP,
    data.frame(node = names(X), X = X, iter = i, estimate = "median")
  )
}

# Summarize bootstrap results
KMSUMMARY <- KMBOOTSTRAP %>% 
  group_by(estimate, node) %>% 
  summarise(bootstrap_mean = mean(X), bootstrap_stddev = sd(X), .groups = "drop")

ACESUMMARY <- ACEBOOTSTRAP %>% 
  group_by(estimate, node) %>% 
  summarise(meanAce = mean(ace), stdev = sd(ace), .groups = "drop")

# Create color palette for visualization
Lab.palette <- colorRampPalette(RColorBrewer::brewer.pal(10, "Spectral"))

# Create visualization of median lifespan
pdf("plots/Figure_median_bootstrap.pdf", width = 16, height = 8)
par(mfrow = c(1, 2), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))

# Plot median lifespan
tipvalues <- KMSUMMARY %>% filter(estimate == "median")
ancvalues <- ACESUMMARY %>% filter(estimate == "median")
X <- tipvalues$bootstrap_mean
names(X) <- tipvalues$node
X <- X[primate_tree_subset$tip.label]
limits <- as.numeric(quantile(X, c(0.05, 0.9)))
cmapout <- contMap(primate_tree_subset, X, lims = limits, plot = FALSE)

plot(setMap(cmapout, colors = Lab.palette(100)), 
     lwd = 7, 
     sig = 1, 
     node.numbers = TRUE, 
     outline = TRUE, 
     font = 1, 
     offset = 1.5, 
     cex = 0.8, 
     ftype = "reg", 
     fsize = 0.9, 
     edge.width = 8,
     leg.txt = "Median Lifespan")

nodelabel <- paste0(round(ancvalues$meanAce, 0))
nodelabels(nodelabel, cex = 1.0, bg = "#FFFFFF", frame = "rect") 
tiplabels(round(X, 0), adj = c(-0.2, 0.5), bg = "#FFFFFF", cex = 0.7, frame = "rect")

# Plot maximum lifespan
tipvalues <- KMSUMMARY %>% filter(estimate == "maximum")
ancvalues <- ACESUMMARY %>% filter(estimate == "maximum")
X <- tipvalues$bootstrap_mean
names(X) <- tipvalues$node
X <- X[primate_tree_subset$tip.label]
limits <- as.numeric(quantile(X, c(0.05, 0.9)))
cmapout <- contMap(primate_tree_subset, X, lims = limits, plot = FALSE)

plot(setMap(cmapout, colors = Lab.palette(100)), 
     lwd = 7, 
     sig = 1, 
     node.numbers = TRUE, 
     outline = TRUE, 
     font = 1, 
     offset = 1.5, 
     cex = 0.8, 
     ftype = "reg", 
     fsize = 0.9, 
     edge.width = 8,
     leg.txt = "Maximum Lifespan")

nodelabel <- paste0(round(ancvalues$meanAce, 0))
nodelabels(nodelabel, cex = 1.0, bg = "#FFFFFF", frame = "rect") 
tiplabels(round(X, 0), adj = c(-0.2, 0.5), bg = "#FFFFFF", cex = 0.7, frame = "rect")
dev.off()

# Create a second visualization focusing on node labels
pdf("plots/Figure_median_bootstrap_2.pdf", width = 16, height = 8)
par(mfrow = c(1, 2), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
plot(setMap(cmapout, colors = Lab.palette(100)), 
     lwd = 0, 
     sig = 1, 
     node.numbers = TRUE, 
     outline = TRUE, 
     font = 1, 
     offset = 1.5, 
     cex = 0.8, 
     ftype = "reg", 
     fsize = 0.9, 
     edge.width = 8)

nodelabel <- ancvalues$node
nodelabels(nodelabel, cex = 1.0, bg = "#FFFFFF", frame = "rect") 
dev.off()
