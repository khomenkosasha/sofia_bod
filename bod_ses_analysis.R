
rm(list=ls()) # clear environment

# load packages
library(tidyverse)
library(sf) 
library(raster)
library(factoextra)
library(psych)
library(car)
library(FSA)
library(reshape2)
library(ggspatial)


#--------------------------------
# SES data - ANOVA
#--------------------------------

### Income by districts
sofia_income <- st_read("Baseline_BOD/income_20170000.geojson")
sofia_income$income_group <- cut(sofia_income$sr_god_dohod, breaks=c(quantile(sofia_income$sr_god_dohod, probs = seq(0, 1, by = 0.20))),labels=c("Very Low","Low","Medium","High", "Very High"), include.lowest = T)

### SES variables from the census
sofia_ses_poly <- st_read("Baseline_BOD/Clean_data/Sofia_SES_poly.geojson")

# Join Income and SES variables
sofia_income <- st_transform(sofia_income, st_crs(sofia_ses_poly))
sofia_ses_poly <- st_join(sofia_income, sofia_ses_poly)

# Remove duplicates based on area coverage
sofia_ses_poly[duplicated(sofia_ses_poly$Poly_ID),]$Poly_ID
sofia_ses_poly$area <- st_area(sofia_ses_poly)
sofia_ses_poly <- sofia_ses_poly %>% group_by(Poly_ID) %>% filter(area==max(area))
sofia_ses_poly <- sofia_ses_poly %>% drop_na(Poly_ID)
sofia_ses_poly$geometry <- NULL
sofia_ses_poly$area <- NULL
summary(sofia_ses_poly)

# Population variables: age, sex 
# sofia_pop <- st_read("Baseline_BOD/Clean_data/sofia_pop_poly.geojson")
# sofia_pop$geometry <- NULL
# # proportion of elderly
# sofia_pop$per_elderly <- ((sofia_pop$female_70. + sofia_pop$male_70.) / sofia_pop$Total) * 100
# summary(sofia_pop$per_elderly)
# sofia_ses_poly <- merge(sofia_ses_poly, sofia_pop[c("Poly_ID", "per_elderly")], by = "Poly_ID")

# Add exposures here
sofia_airp <- st_read("Baseline_BOD/Clean_data/sofia_airp_poly.geojson")
sofia_airp$geometry <- NULL
sofia_noise <- st_read("Baseline_BOD/Clean_data/sofia_noise_poly.geojson")
sofia_noise$geometry <- NULL
sofia_gs <- st_read("Baseline_BOD/Clean_data/sofia_green_poly.geojson")
sofia_gs$geometry <- NULL
sofia_uhi <- st_read("Baseline_BOD/Clean_data/sofia_uhi_poly.geojson")
sofia_uhi$geometry <- NULL
# Merge with SES data
sofia_ses_poly <- merge(sofia_ses_poly, sofia_airp[c("Poly_ID", "pm25", "no2", "o3")], by = "Poly_ID")
sofia_ses_poly <- merge(sofia_ses_poly, sofia_noise[c("Poly_ID", "noise_popw")], by = "Poly_ID")
sofia_ses_poly <- merge(sofia_ses_poly, sofia_gs[c("Poly_ID", "ndvi", "per_green")], by = "Poly_ID")
sofia_ses_poly <- merge(sofia_ses_poly, sofia_uhi[c("Poly_ID", "uhi_summer")], by = "Poly_ID")


### Factor analysis PCA to understand how the variables are interrelated
# Create matrix
ses_matrix <- sofia_ses_poly[c(6, 15:28)]
names(ses_matrix) <- c("Annual Household Income", "% Higher Education", "% Secondary Education", "% Basic Education", "% Primary and lower Education", "% Children", "% Employed", "% Unemployed", "PM2.5", "NO2", "O3", "Noise", "NDVI", "%GA", "UHI")
cor_matrix <- cor(ses_matrix)
#cor_matrix[abs(cor_matrix) < 0.5] <- NA 
print(cor_matrix)

# Correlation plot
cor_melted <- melt(cor_matrix)

# Create the heatmap plot
ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#008B8B", high = "coral", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name = "Correlation") +
  theme_bw() +
  coord_fixed() +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = "", y = "", title = "Correlation SES and environmental variables")
ggsave("Baseline_BOD/Correlation_SES_env.png")

# Factor analysis
ses_matrix <- sofia_ses_poly[c(6, 15:21)]
names(ses_matrix) <- c("Annual Household Income", "% Higher Education", "% Secondary Education", "% Basic Education", "% Primary and lower Education", "% Children", "% Employed", "% Unemployed")
ses_scaled <- scale(ses_matrix) # scale data
fa_result <- fa(ses_scaled, nfactors = 3, rotate = "varimax")
print(fa_result)
fa_result$loadings

# PCA
res.pca <- prcomp(ses_matrix, scale = TRUE)
fviz_eig(res.pca)

# Plot individuals
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             geom = "point",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)     # Avoid text overlapping

# Plot contributions
contrib_dim1 <- get_pca_var(res.pca)$contrib[, 1]
fviz_pca_var(res.pca,
             col.var = contrib_dim1, # Color by contributions to the PC1
             gradient.cols = c("#0066CC", "#FFCC00", "#FF5733"),  # Refined color gradient
             repel = TRUE, 
             legend.title = "Contribution", # Improve legend title
             title = "PCA SES variables contributions", # Title of the plot
             axes.lab.size = 12,   # Adjust axis label size
             title.size = 14,      # Adjust title size
             label.size = 4,       # Adjust variable label size
             axis.textsize = 10)
ggsave("Baseline_BOD/PCA_SES_variables.png")

# PCA results
# Eigenvalues
eig.val <- get_eigenvalue(res.pca)
eig.val

# Results for Variables
res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
# Results for individuals
res.ind <- get_pca_ind(res.pca)
pca_coord <- res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

# Assign coordinates to each neighborhood (Dimension 1)
# sofia_ses_poly$ses_index <- as.vector(pca_coord[, 1])

# Create SES index
ses_matrix <- sofia_ses_poly[c(15:17, 21)]
ses_matrix$H17_1_sum <- 100 - ses_matrix$H17_1_sum
ses_scaled <- scale(ses_matrix) # z-score normalization
sofia_ses_poly$ses_index <- rowMeans(ses_scaled)
# weighted mean based on PCA contributions
#sofia_ses_poly$ses_index <- (ses_scaled[, 1] * 25.690858 + ses_scaled[, 2] * 16.487471 + ses_scaled[, 3] * 11.271821 + ses_scaled[, 4] * 15.050757) / (25.690858 + 16.487471 + 11.271821 + 15.050757)
plot(density(sofia_ses_poly$ses_index))
summary(sofia_ses_poly$ses_index)
# centered around 0 with higher values indicating worse SES


#--------------------------------
# Assign SES quartiles
# Using the following variables: higher edu (proxy for income, employment), secondary edu (proxy for unemployment), and basic education (proxy low edu, children)
sofia_ses_poly$ses_group <- cut(sofia_ses_poly$ses_index, breaks=c(quantile(sofia_ses_poly$ses_index, probs = seq(0, 1, by = 0.25))),labels=c("High SES","Medium-High SES","Medium-Low SES","Low SES"), include.lowest = T)
table(sofia_ses_poly$ses_group)
write.csv(sofia_ses_poly, "Baseline_BOD/Clean_data/Sofia_SESindex_poly.csv", quote = T)


#--------------------------------
### PM2.5 and NO2 levels
sofia_airp <- st_read("Baseline_BOD/Clean_data/sofia_airp_poly.geojson")
sofia_income_airp <- merge(sofia_airp, sofia_ses_poly[c("Poly_ID", "ses_group")], by = "Poly_ID")

# Plot SES index
ggplot(data = sofia_income_airp) +
  geom_sf(aes(fill = ses_group), size = 0.1) + 
  scale_fill_brewer(palette = "BuPu", name = "SES quartiles", guide = guide_legend(nrow = 2, keywidth = 1.5, keyheight = 1.5)) + 
  theme_bw() + # Use a minimal theme
  labs(title = "SES index") +
  theme(legend.position = "bottom", legend.title = element_text(size = 12), legend.text = element_text(size = 12), 
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        panel.grid = element_blank()) + 
  annotation_scale(location = "bl", width_hint = 0.2) +  # Add scale bar at bottom left
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering())

ggsave("Baseline_BOD/SES_index.png")

# Plot by SES group
theme_set(theme_bw())
p1 <- ggplot(sofia_income_airp, aes(x=ses_group, y=pm25)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 20.3, color = "red", linetype = "dashed") +
  labs(title =  expression("PM"[2.5]*" levels by SES group"), x = "", y = expression(paste(mu, "g/m"^3)))
p2 <- ggplot(sofia_income_airp, aes(x=ses_group, y=no2)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 25.8, color = "red", linetype = "dashed") + 
  labs(title = expression("NO"[2]*" levels by SES group"), x = "", y = expression(paste(mu, "g/m"^3)))
p3 <- ggplot(sofia_income_airp, aes(x=ses_group, y=o3)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 62.6, color = "red", linetype = "dashed") + 
  labs(title = expression("O"[3]*" levels by SES group"), x = "", y = expression(paste(mu, "g/m"^3)))

#alternative NO2 data
ggplot(sofia_income_airp, aes(x=ses_group, y=no2_zz)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 31.9, color = "red", linetype = "dashed") + 
  labs(title = expression("NO"[2]*" levels by SES group"), x = "SES group", y = expression(paste(mu, "g/m"^3)))

# Exposure levels by SES group
tapply(sofia_income_airp$pm25, sofia_income_airp$ses_group, summary)
tapply(sofia_income_airp$no2, sofia_income_airp$ses_group, summary)
tapply(sofia_income_airp$no2_zz, sofia_income_airp$ses_group, summary)
tapply(sofia_income_airp$o3, sofia_income_airp$ses_group, summary)

# ANOVA and post-hoc Tukey test
# PM2.5
data.av <- aov(sofia_income_airp$pm25 ~ sofia_income_airp$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_airp$pm25))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_airp$pm25 ~ sofia_income_airp$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_airp$pm25 ~ sofia_income_airp$ses_group)
dunnTest(sofia_income_airp$pm25 ~ sofia_income_airp$ses_group, method="bonferroni")
#pairwise.wilcox.test(sofia_income_airp$pm25, sofia_income_airp$income_group,p.adjust.method = "BH")

# NO2
data.av <- aov(sofia_income_airp$no2 ~ sofia_income_airp$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_airp$no2))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_airp$no2 ~ sofia_income_airp$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_airp$no2 ~ sofia_income_airp$ses_group)
dunnTest(sofia_income_airp$no2 ~ sofia_income_airp$ses_group, method="bonferroni")
#pairwise.wilcox.test(sofia_income_airp$no2, sofia_income_airp$income_group,p.adjust.method = "BH")

# O3
data.av <- aov(sofia_income_airp$o3 ~ sofia_income_airp$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_airp$o3))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_airp$o3 ~ sofia_income_airp$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_airp$o3 ~ sofia_income_airp$ses_group)
dunnTest(sofia_income_airp$o3 ~ sofia_income_airp$ses_group, method="bonferroni")


#--------------------------------
### Environmental noise
sofia_noise <- st_read("Baseline_BOD/Clean_data/sofia_noise_poly.geojson")
sofia_income_noise <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_noise, by = "Poly_ID")

# Plot by SES group
theme_set(theme_bw())
p4 <- ggplot(sofia_income_noise, aes(x=ses_group, y=noise_popw)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 62.2, color = "red", linetype = "dashed") +
  labs(title = "Road-traffic noise levels by SES group", x = "", y = expression(paste("dB(A) Lday")))

# Exposure levels by SES group
tapply(sofia_income_noise$noise_popw, sofia_income_noise$ses_group, summary)

# ANOVA and post-hoc Tukey test
# PM2.5
data.av <- aov(sofia_income_noise$noise_popw ~ sofia_income_noise$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_noise$noise_popw))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_noise$noise_popw ~ sofia_income_noise$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_noise$noise_popw ~ sofia_income_noise$ses_group)
dunnTest(sofia_income_noise$noise_popw ~ sofia_income_noise$ses_group, method="bonferroni")
#pairwise.wilcox.test(sofia_income_airp$pm25, sofia_income_airp$income_group,p.adjust.method = "BH")


#--------------------------------
### Green space levels
sofia_green <- st_read("Baseline_BOD/Clean_data/sofia_green_poly.geojson")
sofia_income_green <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_green, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
p5 <- ggplot(sofia_income_green, aes(x=ses_group, y=ndvi)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 0.397, color = "red", linetype = "dashed") +
  labs(title = "NDVI levels by SES group", x = "", y = "NDVI")

p6 <- ggplot(sofia_income_green, aes(x=ses_group, y=per_green)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 14.4, color = "red", linetype = "dashed") +
  labs(title = "% green area (GA) by SES group", x = "", y = "% GA")

# Exposure levels by income group
tapply(sofia_income_green$ndvi, sofia_income_green$ses_group, summary)
tapply(sofia_income_green$per_green, sofia_income_green$ses_group, summary)

# ANOVA and post-hoc Tukey test
# NDVI
data.av <- aov(sofia_income_green$ndvi ~ sofia_income_green$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_green$ndvi))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_green$ndvi ~ sofia_income_green$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_green$ndvi ~ sofia_income_green$ses_group)
dunnTest(sofia_income_green$ndvi ~ sofia_income_green$ses_group, method="bonferroni")
#pairwise.wilcox.test(sofia_income_green$ndvi, sofia_income_green$income_group,p.adjust.method = "BH")

# % green space
data.av <- aov(sofia_income_green$per_green ~ sofia_income_green$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_green$per_green))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_green$per_green ~ sofia_income_green$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_green$per_green ~ sofia_income_green$ses_group)
dunnTest(sofia_income_green$per_green ~ sofia_income_green$ses_group, method="bonferroni")
#pairwise.wilcox.test(sofia_income_green$ndvi, sofia_income_green$income_group,p.adjust.method = "BH")


#--------------------------------
### UHI effect
sofia_uhi <- st_read("Baseline_BOD/Clean_data/sofia_uhi_poly.geojson")
sofia_income_uhi <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_uhi, by = "Poly_ID")

# Plot by income group
p7 <- ggplot(sofia_income_uhi, aes(x=ses_group, y=uhi_summer)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 2.5, color = "red", linetype = "dashed") +
  labs(title = "UHI effect by SES group", x = "", y = "UHI (ºC)")

# Exposure levels by income group
tapply(sofia_income_uhi$uhi_summer, sofia_income_uhi$ses_group, summary)

# ANOVA and post-hoc Tukey test
# UHI
data.av <- aov(sofia_income_uhi$uhi_summer ~ sofia_income_uhi$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_uhi$uhi_summer))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_uhi$uhi_summer ~ sofia_income_uhi$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_uhi$uhi_summer ~ sofia_income_uhi$ses_group)
dunnTest(sofia_income_uhi$uhi_summer ~ sofia_income_uhi$ses_group, method="bonferroni")


### Export plot with results
library(ggpubr)
ggarrange(p1, p2, p4, p5, p6, p7, nrow = 3, ncol = 2)
ggsave("Baseline_BOD/Sofia_impacts_SES_exposures.png")


#--------------------------------
### PM2.5 mortality impacts
sofia_airp <- st_read("Baseline_BOD/Results/res_pm25_nat_mort.geojson")
sofia_income_airp <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_airp, by = "Poly_ID")
#quantile(sofia_income_airp$att_rate, 0.99)
#sofia_income_airp[sofia_income_airp$att_rate > 456.7762, ]$att_rate <- NA

# Plot by income group
theme_set(theme_bw())
p1 <- ggplot(sofia_income_airp, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 217.4, color = "red", linetype = "dashed") + 
  labs(title = expression("PM"[2.5]*" attributable mortality by SES group"), x = "", y = "Rate per 100,000")

# Mortality rate by income group
tapply(sofia_income_airp$att_rate, sofia_income_airp$ses_group, summary)

# ANOVA and post-hoc Tukey
# PM2.5
data.av <- aov(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_airp$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
dunnTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group, method="bonferroni")
#pairwise.wilcox.test(sofia_income_airp$att_rate, sofia_income_airp$income_group,p.adjust.method = "BH")


#######################
### PM2.5 Hypertension
sofia_airp <- st_read("Baseline_BOD/Results/res_pm25_hypertension.geojson")
sofia_income_airp <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_airp, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
ggplot(sofia_income_airp, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 861, color = "red", linetype = "dashed") + 
  labs(title = expression("PM"[2.5]*" Hypertension incidence by SES group"), x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_airp$att_rate, sofia_income_airp$ses_group, summary)

# ANOVA and post-hoc Tukey
# PM2.5
data.av <- aov(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_airp$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
dunnTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group, method="bonferroni")

### PM2.5 Hypertension SES adjusted
sofia_airp <- st_read("Baseline_BOD/Results/res_pm25_SES_hypertension.geojson")
sofia_income_airp <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_airp, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
ggplot(sofia_income_airp, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 882.6, color = "red", linetype = "dashed") + 
  labs(title = expression("PM"[2.5]*" Hypertension incidence by SES group"), x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_airp$att_rate, sofia_income_airp$ses_group, summary)

# ANOVA and post-hoc Tukey
# PM2.5
data.av <- aov(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_airp$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
dunnTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group, method="bonferroni")


#######################
### PM2.5 IHD
sofia_airp <- st_read("Baseline_BOD/Results/res_pm25_ihd.geojson")
sofia_income_airp <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_airp, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
p2 <- ggplot(sofia_income_airp, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 232, color = "red", linetype = "dashed") + 
  labs(title = expression("PM"[2.5]*" IHD incidence by SES group"), x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_airp$att_rate, sofia_income_airp$ses_group, summary)

# ANOVA and post-hoc Tukey
# PM2.5
data.av <- aov(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_airp$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
dunnTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group, method="bonferroni")

### PM2.5 IHD SES adjusted
sofia_airp <- st_read("Baseline_BOD/Results/res_pm25_SES_ihd.geojson")
sofia_income_airp <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_airp, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
p2 <- ggplot(sofia_income_airp, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 239.3, color = "red", linetype = "dashed") + 
  labs(title = expression("PM"[2.5]*" IHD incidence by SES group"), x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_airp$att_rate, sofia_income_airp$ses_group, summary)

# ANOVA and post-hoc Tukey
# PM2.5
data.av <- aov(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_airp$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
dunnTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group, method="bonferroni")


#######################
### PM2.5 stroke
sofia_airp <- st_read("Baseline_BOD/Results/res_pm25_stroke.geojson")
sofia_income_airp <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_airp, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
ggplot(sofia_income_airp, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 148, color = "red", linetype = "dashed") + 
  labs(title = expression("PM"[2.5]*" Stroke incidence by SES group"), x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_airp$att_rate, sofia_income_airp$ses_group, summary)

# ANOVA and post-hoc Tukey
# PM2.5
data.av <- aov(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_airp$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
dunnTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group, method="bonferroni")

### PM2.5 Stroke SES adjusted
sofia_airp <- st_read("Baseline_BOD/Results/res_pm25_SES_stroke.geojson")
sofia_income_airp <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_airp, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
ggplot(sofia_income_airp, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 151, color = "red", linetype = "dashed") + 
  labs(title = expression("PM"[2.5]*" Stroke incidence by SES group"), x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_airp$att_rate, sofia_income_airp$ses_group, summary)

# ANOVA and post-hoc Tukey
# PM2.5
data.av <- aov(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_airp$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
dunnTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group, method="bonferroni")


#######################
### PM2.5 asthma in children
sofia_airp <- st_read("Baseline_BOD/Results/res_pm25_asthma_children.geojson")
sofia_income_airp <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_airp, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
ggplot(sofia_income_airp, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 207, color = "red", linetype = "dashed") + 
  labs(title = expression("PM"[2.5]*" asthma incidence in children by SES group"), x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_airp$att_rate, sofia_income_airp$ses_group, summary)

# ANOVA and post-hoc Tukey
# PM2.5
data.av <- aov(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_airp$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
dunnTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group, method="bonferroni")


#######################
### PM2.5 COPD
sofia_airp <- st_read("Baseline_BOD/Results/res_pm25_copd.geojson")
sofia_income_airp <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_airp, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
p3 <- ggplot(sofia_income_airp, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 149, color = "red", linetype = "dashed") + 
  labs(title = expression("PM"[2.5]*" COPD incidence by SES group"), x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_airp$att_rate, sofia_income_airp$ses_group, summary)

# ANOVA and post-hoc Tukey
# PM2.5
data.av <- aov(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_airp$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
dunnTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group, method="bonferroni")

### PM2.5 COPD SES adjusted
sofia_airp <- st_read("Baseline_BOD/Results/res_pm25_SES_copd.geojson")
sofia_income_airp <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_airp, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
p3 <- ggplot(sofia_income_airp, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 151.9, color = "red", linetype = "dashed") + 
  labs(title = expression("PM"[2.5]*" COPD incidence by SES group"), x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_airp$att_rate, sofia_income_airp$ses_group, summary)

# ANOVA and post-hoc Tukey
# PM2.5
data.av <- aov(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_airp$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
dunnTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group, method="bonferroni")


#######################
### PM2.5 diabetes
sofia_airp <- st_read("Baseline_BOD/Results/res_pm25_diabetes.geojson")
sofia_income_airp <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_airp, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
ggplot(sofia_income_airp, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 140, color = "red", linetype = "dashed") + 
  labs(title = expression("PM"[2.5]*" Diabetes (type II) incidence by SES group"), x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_airp$att_rate, sofia_income_airp$ses_group, summary)

# ANOVA and post-hoc Tukey
# PM2.5
data.av <- aov(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_airp$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
dunnTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group, method="bonferroni")

### PM2.5 diabetes SES adjusted
sofia_airp <- st_read("Baseline_BOD/Results/res_pm25_SES_diabetes.geojson")
sofia_income_airp <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_airp, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
ggplot(sofia_income_airp, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 142.7, color = "red", linetype = "dashed") + 
  labs(title = expression("PM"[2.5]*" Diabetes (type II) incidence by SES group"), x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_airp$att_rate, sofia_income_airp$ses_group, summary)

# ANOVA and post-hoc Tukey
# PM2.5
data.av <- aov(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_airp$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
dunnTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group, method="bonferroni")


#--------------------------------
### NO2 mortality impacts
sofia_airp <- st_read("Baseline_BOD/Results/res_no2_nat_mort.geojson")
sofia_income_airp <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_airp, by = "Poly_ID")
#quantile(sofia_income_airp$att_rate, 0.99)
#sofia_income_airp[sofia_income_airp$att_rate > 294.9395, ]$att_rate <- NA

# Plot by income group
theme_set(theme_bw())
p4 <- ggplot(sofia_income_airp, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 131.4, color = "red", linetype = "dashed") +
  labs(title = expression("NO"[2]*" attributable mortality by SES group"), x = "", y = "Rate per 100,000")

# Mortality rate by income group
tapply(sofia_income_airp$att_rate, sofia_income_airp$ses_group, summary)

# ANOVA and post-hoc Tukey
# NO2
data.av <- aov(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_airp$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
dunnTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group, method="bonferroni")
#pairwise.wilcox.test(sofia_income_airp$att_rate, sofia_income_airp$income_group,p.adjust.method = "BH")


###########################
### NO2 asthma in adults
sofia_airp <- st_read("Baseline_BOD/Results/res_no2_asthma_adults.geojson")
sofia_income_airp <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_airp, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
p5 <- ggplot(sofia_income_airp, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 23, color = "red", linetype = "dashed") +
  labs(title = expression("NO"[2]*" asthma incidence in adults by SES group"), x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_airp$att_rate, sofia_income_airp$ses_group, summary)

# ANOVA and post-hoc Tukey
# PM2.5
data.av <- aov(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_airp$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
dunnTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group, method="bonferroni")

### NO2 asthma in adults SES adjusted
sofia_airp <- st_read("Baseline_BOD/Results/res_no2_SES_asthma_adults.geojson")
sofia_income_airp <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_airp, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
p5 <- ggplot(sofia_income_airp, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 24.5, color = "red", linetype = "dashed") +
  labs(title = expression("NO"[2]*" asthma incidence in adults by SES group"), x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_airp$att_rate, sofia_income_airp$ses_group, summary)

# ANOVA and post-hoc Tukey
# PM2.5
data.av <- aov(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_airp$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
dunnTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group, method="bonferroni")


###########################
### NO2 asthma in children
sofia_airp <- st_read("Baseline_BOD/Results/res_no2_asthma_children.geojson")
sofia_income_airp <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_airp, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
ggplot(sofia_income_airp, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 84, color = "red", linetype = "dashed") +
  labs(title = expression("NO"[2]*" asthma incidence in children by SES group"), x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_airp$att_rate, sofia_income_airp$ses_group, summary)

# ANOVA and post-hoc Tukey
# PM2.5
data.av <- aov(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_airp$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
dunnTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group, method="bonferroni")


###########################
### NO2 ALRI in children 
sofia_airp <- st_read("Baseline_BOD/Results/res_no2_alri_children.geojson")
sofia_income_airp <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_airp, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
ggplot(sofia_income_airp, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 514, color = "red", linetype = "dashed") +
  labs(title = expression("NO"[2]*" ALRI incidence in children by SES group"), x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_airp$att_rate, sofia_income_airp$ses_group, summary)

# ANOVA and post-hoc Tukey
# PM2.5
data.av <- aov(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_airp$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
dunnTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group, method="bonferroni")


#--------------------------------
### PM2.5 and NO2 multiexposure mortality impacts
sofia_airp <- st_read("Baseline_BOD/Results/res_pm25_nat_mort_multiexp.geojson")
sofia_income_airp <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_airp, by = "Poly_ID")
#quantile(sofia_income_airp$att_rate, 0.99)
#sofia_income_airp[sofia_income_airp$att_rate > 503.4286, ]$att_rate <- NA

# Plot by income group
theme_set(theme_bw())
p6 <- ggplot(sofia_income_airp, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 244.1, color = "red", linetype = "dashed") +
  labs(title = expression("PM"[2.5]*" and NO"[2]*" attributable mortality by SES group"), x = "", y = "Rate per 100,000")

# Mortality rate by income group
tapply(sofia_income_airp$att_rate, sofia_income_airp$ses_group, summary)

# ANOVA and post-hoc Tukey
# Multi-exposure
data.av <- aov(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_airp$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
dunnTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group, method="bonferroni")
#pairwise.wilcox.test(sofia_income_airp$att_rate, sofia_income_airp$income_group,p.adjust.method = "BH")


#--------------------------------
### O3 mortality impacts
sofia_airp <- st_read("Baseline_BOD/Results/res_o3_resp_mort.geojson")
sofia_income_airp <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_airp, by = "Poly_ID")
#quantile(sofia_income_airp$att_rate, 0.99)
#sofia_income_airp[sofia_income_airp$att_rate > 294.9395, ]$att_rate <- NA

# Plot by income group
theme_set(theme_bw())
ggplot(sofia_income_airp, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 4.9, color = "red", linetype = "dashed") +
  labs(title = expression("O"[3]*" attributable mortality by SES group"), x = "", y = "Rate per 100,000")

# Mortality rate by income group
tapply(sofia_income_airp$att_rate, sofia_income_airp$ses_group, summary)

# ANOVA and post-hoc Tukey
# NO2
data.av <- aov(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_airp$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group)
dunnTest(sofia_income_airp$att_rate ~ sofia_income_airp$ses_group, method="bonferroni")
#pairwise.wilcox.test(sofia_income_airp$att_rate, sofia_income_airp$income_group,p.adjust.method = "BH")


#--------------------------------
### Noise mortality impacts
sofia_noise <- st_read("Baseline_BOD/Results/res_noise_nat_mort.geojson")
sofia_income_noise <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_noise, by = "Poly_ID")
#quantile(sofia_income_noise$att_rate, 0.99)
#sofia_income_noise[sofia_income_noise$att_rate > 199.7877, ]$att_rate <- NA

# Plot by income group
theme_set(theme_bw())
p7 <- ggplot(sofia_income_noise, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 101, color = "red", linetype = "dashed") +
  labs(title = "Noise attributable mortality by SES group", x = "", y = "Rate per 100,000")

# Mortality rate by income group
tapply(sofia_income_noise$att_rate, sofia_income_noise$ses_group, summary)

# ANOVA and post-hoc Tukey
# NO2
data.av <- aov(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_noise$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group)
dunnTest(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group, method="bonferroni")
#pairwise.wilcox.test(sofia_income_airp$att_rate, sofia_income_airp$income_group,p.adjust.method = "BH")


##########################
### Noise high annoyance 
sofia_noise <- st_read("Baseline_BOD/Results/res_noise_annoyance.geojson")
sofia_income_noise <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_noise, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
ggplot(sofia_income_noise, aes(x=ses_group, y=mean)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 22.1, color = "red", linetype = "dashed") +
  labs(title = "High noise annoyance (%) by SES group", x = "", y = "% of population")

# Incidence rate by income group
tapply(sofia_income_noise$mean, sofia_income_noise$ses_group, summary)

# ANOVA and post-hoc Tukey
# NO2
data.av <- aov(sofia_income_noise$mean ~ sofia_income_noise$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_noise$mean, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_noise$mean ~ sofia_income_noise$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_noise$mean ~ sofia_income_noise$ses_group)
dunnTest(sofia_income_noise$mean ~ sofia_income_noise$ses_group, method="bonferroni")


##########################
### Noise high sleep disturbance 
sofia_noise <- st_read("Baseline_BOD/Results/res_noise_sleep_dist.geojson")
sofia_income_noise <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_noise, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
ggplot(sofia_income_noise, aes(x=ses_group, y=mean)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 7.5, color = "red", linetype = "dashed") +
  labs(title = "High sleep disturbance (%) by SES group", x = "", y = "% of population")

# Incidence rate by income group
tapply(sofia_income_noise$mean, sofia_income_noise$ses_group, summary)

# ANOVA and post-hoc Tukey
# NO2
data.av <- aov(sofia_income_noise$mean ~ sofia_income_noise$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_noise$mean, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_noise$mean ~ sofia_income_noise$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_noise$mean ~ sofia_income_noise$ses_group)
dunnTest(sofia_income_noise$mean ~ sofia_income_noise$ses_group, method="bonferroni")


##########################
### Noise IHD impacts 
sofia_noise <- st_read("Baseline_BOD/Results/res_noise_ihd.geojson")
sofia_income_noise <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_noise, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
p8 <- ggplot(sofia_income_noise, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 63, color = "red", linetype = "dashed") +
  labs(title = "Noise IHD incidence by SES group", x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_noise$att_rate, sofia_income_noise$ses_group, summary)

# ANOVA and post-hoc Tukey
# Noise
data.av <- aov(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_noise$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group)
dunnTest(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group, method="bonferroni")

### Noise IHD impacts SES adjusted 
sofia_noise <- st_read("Baseline_BOD/Results/res_noise_SES_ihd.geojson")
sofia_income_noise <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_noise, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
p8 <- ggplot(sofia_income_noise, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 62.5, color = "red", linetype = "dashed") +
  labs(title = "Noise IHD incidence by SES group", x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_noise$att_rate, sofia_income_noise$ses_group, summary)

# ANOVA and post-hoc Tukey
# NO2
data.av <- aov(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_noise$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group)
dunnTest(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group, method="bonferroni")


##########################
### Noise stroke impacts 
sofia_noise <- st_read("Baseline_BOD/Results/res_noise_stroke.geojson")
sofia_income_noise <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_noise, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
ggplot(sofia_income_noise, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 37, color = "red", linetype = "dashed") +
  labs(title = "Noise Stroke incidence by SES group", x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_noise$att_rate, sofia_income_noise$ses_group, summary)

# ANOVA and post-hoc Tukey
# Noise
data.av <- aov(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_noise$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group)
dunnTest(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group, method="bonferroni")

### Noise stroke impacts SES adjusted 
sofia_noise <- st_read("Baseline_BOD/Results/res_noise_SES_stroke.geojson")
sofia_income_noise <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_noise, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
ggplot(sofia_income_noise, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 36.5, color = "red", linetype = "dashed") +
  labs(title = "Noise Stroke incidence by SES group", x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_noise$att_rate, sofia_income_noise$ses_group, summary)

# ANOVA and post-hoc Tukey
# NO2
data.av <- aov(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_noise$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group)
dunnTest(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group, method="bonferroni")


##########################
### Noise diabetes impacts 
sofia_noise <- st_read("Baseline_BOD/Results/res_noise_diabetes.geojson")
sofia_income_noise <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_noise, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
p9 <- ggplot(sofia_income_noise, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 70, color = "red", linetype = "dashed") +
  labs(title = "Noise diabetes (type II) incidence by SES group", x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_noise$att_rate, sofia_income_noise$ses_group, summary)

# ANOVA and post-hoc Tukey
# Noise
data.av <- aov(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_noise$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group)
dunnTest(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group, method="bonferroni")

### Noise diabetes impacts SES adjusted 
sofia_noise <- st_read("Baseline_BOD/Results/res_noise_SES_diabetes.geojson")
sofia_income_noise <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_noise, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
p9 <- ggplot(sofia_income_noise, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 68.7, color = "red", linetype = "dashed") +
  labs(title = "Noise diabetes (type II) incidence by SES group", x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_noise$att_rate, sofia_income_noise$ses_group, summary)

# ANOVA and post-hoc Tukey
# NO2
data.av <- aov(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_noise$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group)
dunnTest(sofia_income_noise$att_rate ~ sofia_income_noise$ses_group, method="bonferroni")


#--------------------------------
### Green space mortality impacts
sofia_green <- st_read("Baseline_BOD/Results/res_ndvi_nat_mort.geojson")
sofia_income_green <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_green, by = "Poly_ID")
#quantile(sofia_income_green$att_rate, 0.99)
#sofia_income_green[sofia_income_green$att_rate > 167.1165, ]$att_rate <- NA

# Plot by income group
theme_set(theme_bw())
p10 <- ggplot(sofia_income_green, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 24, color = "red", linetype = "dashed") +
  labs(title = "Lack of green attributable mortality by SES group", x = "", y = "Rate per 100,000")
ggsave("Baseline_BOD/Sofia_impacts_SES_NDVI_mortality.png")

# Mortality rate by income group
tapply(sofia_income_green$att_rate, sofia_income_green$ses_group, summary)

# ANOVA and post-hoc Tukey test
# Green
data.av <- aov(sofia_income_green$att_rate ~ sofia_income_green$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_green$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_green$att_rate ~ sofia_income_green$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_green$att_rate ~ sofia_income_green$ses_group)
dunnTest(sofia_income_green$att_rate ~ sofia_income_green$ses_group, method="bonferroni")
#pairwise.wilcox.test(sofia_income_green$att_rate, sofia_income_green$income_group,p.adjust.method = "BH")


##############################
### Green hypertension impacts
sofia_green <- st_read("Baseline_BOD/Results/res_ndvi_hypertension.geojson")
sofia_income_green <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_green, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
ggplot(sofia_income_green, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 29, color = "red", linetype = "dashed") +
  labs(title = "Lack of green hypertension incidence by SES group", x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_green$att_rate, sofia_income_green$ses_group, summary)

# ANOVA and post-hoc Tukey test
# Green
data.av <- aov(sofia_income_green$att_rate ~ sofia_income_green$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_green$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_green$att_rate ~ sofia_income_green$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_green$att_rate ~ sofia_income_green$ses_group)
dunnTest(sofia_income_green$att_rate ~ sofia_income_green$ses_group, method="bonferroni")

### Green hypertension impacts SES adjusted
sofia_green <- st_read("Baseline_BOD/Results/res_ndvi_SES_hypertension.geojson")
sofia_income_green <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_green, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
ggplot(sofia_income_green, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 28.8, color = "red", linetype = "dashed") +
  labs(title = "Lack of green hypertension incidence by SES group", x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_green$att_rate, sofia_income_green$ses_group, summary)

# ANOVA and post-hoc Tukey test
# Green
data.av <- aov(sofia_income_green$att_rate ~ sofia_income_green$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_green$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_green$att_rate ~ sofia_income_green$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_green$att_rate ~ sofia_income_green$ses_group)
dunnTest(sofia_income_green$att_rate ~ sofia_income_green$ses_group, method="bonferroni")


##############################
### Green stroke impacts
sofia_green <- st_read("Baseline_BOD/Results/res_ndvi_stroke.geojson")
sofia_income_green <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_green, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
p11 <- ggplot(sofia_income_green, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 5, color = "red", linetype = "dashed") +
  labs(title = "Lack of green stroke incidence by SES group", x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_green$att_rate, sofia_income_green$ses_group, summary)

# ANOVA and post-hoc Tukey test
# Green
data.av <- aov(sofia_income_green$att_rate ~ sofia_income_green$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_green$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_green$att_rate ~ sofia_income_green$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_green$att_rate ~ sofia_income_green$ses_group)
dunnTest(sofia_income_green$att_rate ~ sofia_income_green$ses_group, method="bonferroni")

### Green stroke impacts SES adjusted
sofia_green <- st_read("Baseline_BOD/Results/res_ndvi_SES_stroke.geojson")
sofia_income_green <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_green, by = "Poly_ID")

# Plot by income group
theme_set(theme_bw())
p11 <- ggplot(sofia_income_green, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 5, color = "red", linetype = "dashed") +
  labs(title = "Lack of green stroke incidence by SES group", x = "", y = "Rate per 100,000")

# Incidence rate by income group
tapply(sofia_income_green$att_rate, sofia_income_green$ses_group, summary)

# ANOVA and post-hoc Tukey test
# Green
data.av <- aov(sofia_income_green$att_rate ~ sofia_income_green$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_green$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_green$att_rate ~ sofia_income_green$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_green$att_rate ~ sofia_income_green$ses_group)
dunnTest(sofia_income_green$att_rate ~ sofia_income_green$ses_group, method="bonferroni")


#--------------------------------
### UHI health impacts
### UHI mortality impacts
sofia_uhi <- st_read("Baseline_BOD/Results/res_uhi_allcause_mort.geojson")
sofia_income_uhi <- merge(sofia_ses_poly[c("Poly_ID", "ses_group")], sofia_uhi, by = "Poly_ID")
#quantile(sofia_income_uhi$att_rate, 0.99)
#sofia_income_uhi[sofia_income_uhi$att_rate > 12.7309, ]$att_rate <- NA

# Plot by income group
theme_set(theme_bw())
p12 <- ggplot(sofia_income_uhi, aes(x=ses_group, y=att_rate)) + 
  geom_boxplot(fill = "gray95") + geom_hline(yintercept = 10, color = "red", linetype = "dashed") +
  labs(title = "UHI attributable mortality by SES group", x = "", y = "Rate per 100,000")

# Mortality rate by income group
tapply(sofia_income_uhi$att_rate, sofia_income_uhi$ses_group, summary)

# ANOVA and post-hoc Tukey test
# Green
data.av <- aov(sofia_income_uhi$att_rate ~ sofia_income_uhi$ses_group)
summary(data.av)

# Check for the normality assumption
plot(density(sofia_income_uhi$att_rate, na.rm = T))
plot(data.av, 1) # check homogeneity of variances
leveneTest(sofia_income_uhi$att_rate ~ sofia_income_uhi$ses_group)
plot(data.av, 2) # check normality
# Kruskal wallis test
kruskal.test(sofia_income_uhi$att_rate ~ sofia_income_uhi$ses_group)
dunnTest(sofia_income_uhi$att_rate ~ sofia_income_uhi$ses_group, method="bonferroni")


### Export plot with results
library(ggpubr)

# Mortality outcomes
ggarrange(p1, p4, p6, p7, p10, p12, nrow = 3, ncol = 2)
ggsave("Baseline_BOD/Sofia_impacts_SES_mortality.png")

# Morbidity outcomes
ggarrange(p2, p3, p5, p8, p9, p11, nrow = 3, ncol = 2)
ggsave("Baseline_BOD/Sofia_impacts_SES_morbidity.png")


#--------------------------------
# SES data - Moran I
#--------------------------------

library(spdep)
library(rgdal)
library(tmap)
library(robustHD)


#-----------------------------------
### PM2.5, NO2 and 03 levels
sofia_airp <- st_read("Baseline_BOD/Clean_data/sofia_airp_poly.geojson")
sofia_income_airp <- merge(sofia_airp, sofia_ses_poly[c("Poly_ID", "ses_index")], by = "Poly_ID")
sofia_income_airp$ses_index <- -sofia_income_airp$ses_index

# Plot data
qtm(sofia_income_airp, "pm25")
qtm(sofia_income_airp, "no2")
qtm(sofia_income_airp, "o3")
qtm(sofia_income_airp, "ses_index")

# Spatial weights matrix
centroids <- st_centroid(st_geometry(sofia_income_airp))
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
isolated_areas <- which(card(wm_d) == 0) # Identify areas with zero neighbors
sofia_income_airp <- sofia_income_airp[-isolated_areas, ] # Remove isolated areas from the data
centroids <- centroids[-isolated_areas, ]
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
# Row standardized weights matrix 
rswm_d <- nb2listw(wm_d, style = "W", zero.policy = TRUE)
summary(rswm_d)

### Analysis for PM2.5
# Bivariate Moran I scatter plot
sofia_income_airp$std_ses <- standardize(sofia_income_airp$ses_index)
sofia_income_airp$std_pm25 <- standardize(sofia_income_airp$pm25)
# Lag variable: 
sofia_income_airp$wsr_pm25 <- lag.listw(rswm_d,sofia_income_airp$std_pm25)

# Plot: 
theme_set(theme_bw())
p1 <- ggplot(data = sofia_income_airp, aes(x=std_ses,y=wsr_pm25)) +
  geom_point(color = "grey", alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty =2) +
  ggtitle(expression("SES and lagged PM"[2.5]*"")) + ylab(expression("Standard lagged PM"[2.5]*"")) + xlab("Standard SES index")

### Global Moran I (bivariate)
bivariate_moran <- moran_bv(sofia_income_airp$std_ses, sofia_income_airp$std_pm25, rswm_d, nsim = 500)
print(bivariate_moran)
bivariate_moran$t0
summary(lm(sofia_income_airp$wsr_pm25 ~ sofia_income_airp$std_ses))

### Local Moran I (bivariate)
localMI <- localmoran_bv(sofia_income_airp$std_ses, sofia_income_airp$std_pm25, rswm_d, nsim = 100)
head(localMI)

# Plot results
sofia_income_airp.localMI <- cbind(sofia_income_airp, localMI)

# Moran I
localMI.map <- tm_shape(sofia_income_airp.localMI) +
  tm_fill(col = "Ibvi",
          style = "pretty",
          title = "local moran statistics")+
  tm_borders(alpha = 0.5)
# p values
pvalue.map <- tm_shape(sofia_income_airp.localMI) +
  tm_fill(col = "Pr.z....E.Ibvi..",
          breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
          palette = "-Blues",
          title = "local Moran's I p-values")+
  tm_borders(alpha = 0.5)

tmap_arrange(localMI.map, pvalue.map, asp = 1, ncol = 2)

### LISA cluster map
# Identifying the LISA clusters
quadrant <- vector(mode = "numeric", length = nrow(localMI))
signif <- 0.05 # significance level

# Assign quadrant values (HH, LL, HL, LH)
quadrant[sofia_income_airp.localMI$std_ses > 0 & sofia_income_airp.localMI$wsr_pm25 > 0] <- "High SES-High PM2.5"      
quadrant[sofia_income_airp.localMI$std_ses < 0 & sofia_income_airp.localMI$wsr_pm25 < 0] <- "Low SES-Low PM2.5"      
quadrant[sofia_income_airp.localMI$std_ses < 0 & sofia_income_airp.localMI$wsr_pm25 > 0] <- "Low SES-High PM2.5"
quadrant[sofia_income_airp.localMI$std_ses > 0 & sofia_income_airp.localMI$wsr_pm25 < 0] <- "High SES-Low PM2.5"
quadrant[localMI[,5]>signif] <- "NS"
sofia_income_airp.localMI$quadrant <- quadrant 

# Merge back to full polygon
sofia_income_airp.localMI$geometry <- NULL
sofia_income_airp <- merge(sofia_airp, sofia_income_airp.localMI[c("Poly_ID", "quadrant")], by = "Poly_ID", all.x = T)
sofia_income_airp <- sofia_income_airp %>% mutate(quadrant = replace_na(quadrant, "NS"))
sofia_income_airp$quadrant <- factor(sofia_income_airp$quadrant, levels = c("Low SES-High PM2.5", "High SES-High PM2.5", "High SES-Low PM2.5", "Low SES-Low PM2.5", "NS"))

# Make the plot
theme_set(theme_bw())
w1 <- ggplot(sofia_income_airp) + geom_sf(aes(fill = quadrant), size = 0.1) + 
          scale_fill_manual(name = "", values = c("High SES-High PM2.5" = "#fdae61",  
                                       "Low SES-Low PM2.5" = "#abd9e9", 
                                       "Low SES-High PM2.5" = "#d7191c", 
                                       "High SES-Low PM2.5" = "#2c7bb6", 
                                       "NS" = "#ffffff"), guide = guide_legend(nrow = 2, keywidth = 1.5, keyheight = 1.5)) +
          ggtitle(expression("PM"[2.5]*" exposure")) + 
          theme(legend.position = "bottom", legend.text = element_text(size = 12), 
                axis.text = element_blank(),  # Remove axis text
                axis.ticks = element_blank(), # Remove axis ticks
                panel.grid = element_blank()) + 
  annotation_scale(location = "bl", width_hint = 0.2) +  # Add scale bar at bottom left
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering())


ggsave("Baseline_BOD/Sofia_spatial_analysis_PM2.5.png")


### Analysis for NO2
sofia_airp <- st_read("Baseline_BOD/Clean_data/sofia_airp_poly.geojson")
sofia_income_airp <- merge(sofia_airp, sofia_ses_poly[c("Poly_ID", "ses_index")], by = "Poly_ID")
sofia_income_airp$ses_index <- -sofia_income_airp$ses_index

# Spatial weights matrix
centroids <- st_centroid(st_geometry(sofia_income_airp))
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
isolated_areas <- which(card(wm_d) == 0) # Identify areas with zero neighbors
sofia_income_airp <- sofia_income_airp[-isolated_areas, ] # Remove isolated areas from the data
centroids <- centroids[-isolated_areas, ]
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
# Row standardized weights matrix 
rswm_d <- nb2listw(wm_d, style = "W", zero.policy = TRUE)
summary(rswm_d)

# Bivariate Moran I scatter plot
sofia_income_airp$std_ses <- standardize(sofia_income_airp$ses_index)
sofia_income_airp$std_no2 <- standardize(sofia_income_airp$no2)
# Lag variable: 
sofia_income_airp$wsr_no2 <- lag.listw(rswm_d,sofia_income_airp$std_no2)

# Plot: 
theme_set(theme_bw())
p2 <- ggplot(data = sofia_income_airp, aes(x=std_ses,y=wsr_no2)) +
  geom_point(color = "grey", alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty =2) +
  ggtitle(expression("SES and lagged NO"[2]*"")) + ylab(expression("Standard lagged NO"[2]*"")) + xlab("Standard SES index")

### Global Moran I (bivariate)
bivariate_moran <- moran_bv(sofia_income_airp$std_ses, sofia_income_airp$std_no2, rswm_d, nsim = 500)
print(bivariate_moran)
bivariate_moran$t0
summary(lm(sofia_income_airp$wsr_no2 ~ sofia_income_airp$std_ses))

### Local Moran I (bivariate)
localMI <- localmoran_bv(sofia_income_airp$std_ses, sofia_income_airp$std_no2, rswm_d, nsim = 100)
head(localMI)

# Plot results
sofia_income_airp.localMI <- cbind(sofia_income_airp, localMI)

# Moran I
localMI.map <- tm_shape(sofia_income_airp.localMI) +
  tm_fill(col = "Ibvi",
          style = "pretty",
          title = "local moran statistics")+
  tm_borders(alpha = 0.5)
# p values
pvalue.map <- tm_shape(sofia_income_airp.localMI) +
  tm_fill(col = "Pr.z....E.Ibvi..",
          breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
          palette = "-Blues",
          title = "local Moran's I p-values")+
  tm_borders(alpha = 0.5)

tmap_arrange(localMI.map, pvalue.map, asp = 1, ncol = 2)

### LISA cluster map
# Identifying the LISA clusters
quadrant <- vector(mode = "numeric", length = nrow(localMI))
signif <- 0.05 # significance level

# Assign quadrant values (HH, LL, HL, LH)
quadrant[sofia_income_airp.localMI$std_ses > 0 & sofia_income_airp.localMI$wsr_no2 > 0] <- "High SES-High NO2"      
quadrant[sofia_income_airp.localMI$std_ses < 0 & sofia_income_airp.localMI$wsr_no2 < 0] <- "Low SES-Low NO2"      
quadrant[sofia_income_airp.localMI$std_ses < 0 & sofia_income_airp.localMI$wsr_no2 > 0] <- "Low SES-High NO2"
quadrant[sofia_income_airp.localMI$std_ses > 0 & sofia_income_airp.localMI$wsr_no2 < 0] <- "High SES-Low NO2"
quadrant[localMI[,5]>signif] <- "NS"
sofia_income_airp.localMI$quadrant <- quadrant 

# Merge back to full polygon
sofia_income_airp.localMI$geometry <- NULL
sofia_income_airp <- merge(sofia_airp, sofia_income_airp.localMI[c("Poly_ID", "quadrant")], by = "Poly_ID", all.x = T)
sofia_income_airp <- sofia_income_airp %>% mutate(quadrant = replace_na(quadrant, "NS"))
sofia_income_airp$quadrant <- factor(sofia_income_airp$quadrant, levels = c("Low SES-High NO2", "High SES-High NO2", "High SES-Low NO2", "Low SES-Low NO2", "NS"))

# Make the plot
theme_set(theme_bw())
w3 <- ggplot(sofia_income_airp) + geom_sf(aes(fill = quadrant), size = 0.1) + 
  scale_fill_manual(name = "", values = c("High SES-High NO2" = "#fdae61",  
                                          "Low SES-Low NO2" = "#abd9e9", 
                                          "Low SES-High NO2" = "#d7191c", 
                                          "High SES-Low NO2" = "#2c7bb6", 
                                          "NS" = "#ffffff"), guide = guide_legend(nrow = 2, keywidth = 1.5, keyheight = 1.5)) +
  ggtitle(expression("NO"[2]*" exposure")) + 
  theme(legend.position = "bottom", legend.text = element_text(size = 12), 
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        panel.grid = element_blank()) + 
  annotation_scale(location = "bl", width_hint = 0.2) +  # Add scale bar at bottom left
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering())

ggsave("Baseline_BOD/Sofia_spatial_analysis_NO2.png")


### Analysis for O3
sofia_airp <- st_read("Baseline_BOD/Clean_data/sofia_airp_poly.geojson")
sofia_income_airp <- merge(sofia_airp, sofia_ses_poly[c("Poly_ID", "ses_index")], by = "Poly_ID")
sofia_income_airp$ses_index <- -sofia_income_airp$ses_index

# Spatial weights matrix
centroids <- st_centroid(st_geometry(sofia_income_airp))
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
isolated_areas <- which(card(wm_d) == 0) # Identify areas with zero neighbors
sofia_income_airp <- sofia_income_airp[-isolated_areas, ] # Remove isolated areas from the data
centroids <- centroids[-isolated_areas, ]
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
# Row standardized weights matrix 
rswm_d <- nb2listw(wm_d, style = "W", zero.policy = TRUE)
summary(rswm_d)

# Bivariate Moran I scatter plot
sofia_income_airp$std_ses <- standardize(sofia_income_airp$ses_index)
sofia_income_airp$std_o3 <- standardize(sofia_income_airp$o3)
# Lag variable: 
sofia_income_airp$wsr_o3 <- lag.listw(rswm_d,sofia_income_airp$std_o3)

# Plot: 
theme_set(theme_bw())
p3 <- ggplot(data = sofia_income_airp, aes(x=std_ses,y=wsr_o3)) +
  geom_point(color = "grey", alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty =2) +
  ggtitle(expression("SES and lagged O"[3]*"")) + ylab(expression("Standard lagged O"[3]*"")) + xlab("Standard SES index")

### Global Moran I (bivariate)
bivariate_moran <- moran_bv(sofia_income_airp$std_ses, sofia_income_airp$std_o3, rswm_d, nsim = 500)
print(bivariate_moran)
bivariate_moran$t0
summary(lm(sofia_income_airp$wsr_o3 ~ sofia_income_airp$std_ses))

### Local Moran I (bivariate)
localMI <- localmoran_bv(sofia_income_airp$std_ses, sofia_income_airp$std_o3, rswm_d, nsim = 100)
head(localMI)

# Plot results
sofia_income_airp.localMI <- cbind(sofia_income_airp, localMI)

# Moran I
localMI.map <- tm_shape(sofia_income_airp.localMI) +
  tm_fill(col = "Ibvi",
          style = "pretty",
          title = "local moran statistics")+
  tm_borders(alpha = 0.5)
# p values
pvalue.map <- tm_shape(sofia_income_airp.localMI) +
  tm_fill(col = "Pr.z....E.Ibvi..",
          breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
          palette = "-Blues",
          title = "local Moran's I p-values")+
  tm_borders(alpha = 0.5)

tmap_arrange(localMI.map, pvalue.map, asp = 1, ncol = 2)

### LISA cluster map
# Identifying the LISA clusters
quadrant <- vector(mode = "numeric", length = nrow(localMI))
signif <- 0.05 # significance level

# Assign quadrant values (HH, LL, HL, LH)
quadrant[sofia_income_airp.localMI$std_ses > 0 & sofia_income_airp.localMI$wsr_o3 > 0] <- "High SES-High O3"      
quadrant[sofia_income_airp.localMI$std_ses < 0 & sofia_income_airp.localMI$wsr_o3 < 0] <- "Low SES-Low O3"      
quadrant[sofia_income_airp.localMI$std_ses < 0 & sofia_income_airp.localMI$wsr_o3 > 0] <- "Low SES-High O3"
quadrant[sofia_income_airp.localMI$std_ses > 0 & sofia_income_airp.localMI$wsr_o3 < 0] <- "High SES-Low O3"
quadrant[localMI[,5]>signif] <- "NS"
sofia_income_airp.localMI$quadrant <- quadrant 

# Merge back to full polygon
sofia_income_airp.localMI$geometry <- NULL
sofia_income_airp <- merge(sofia_airp, sofia_income_airp.localMI[c("Poly_ID", "quadrant")], by = "Poly_ID", all.x = T)
sofia_income_airp <- sofia_income_airp %>% mutate(quadrant = replace_na(quadrant, "NS"))
sofia_income_airp$quadrant <- factor(sofia_income_airp$quadrant, levels = c("Low SES-High O3", "High SES-High O3", "High SES-Low O3", "Low SES-Low O3", "NS"))

# Make the plot
theme_set(theme_bw())
w4 <- ggplot(sofia_income_airp) + geom_sf(aes(fill = quadrant), size = 0.1) + 
  scale_fill_manual(name = "", values = c("High SES-High O3" = "#fdae61",  
                                          "Low SES-Low O3" = "#abd9e9", 
                                          "Low SES-High O3" = "#d7191c", 
                                          "High SES-Low O3" = "#2c7bb6", 
                                          "NS" = "#ffffff"), guide = guide_legend(nrow = 2, keywidth = 1.5, keyheight = 1.5)) +
  ggtitle(expression("O"[3]*" exposure")) + 
  theme(legend.position = "bottom", legend.text = element_text(size = 12), 
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        panel.grid = element_blank()) + 
  annotation_scale(location = "bl", width_hint = 0.2) +  # Add scale bar at bottom left
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering())

ggsave("Baseline_BOD/Sofia_spatial_analysis_O3.png")


#-----------------------------------
### Noise levels
sofia_noise <- st_read("Baseline_BOD/Clean_data/sofia_noise_poly.geojson")
sofia_income_noise <- merge(sofia_noise, sofia_ses_poly[c("Poly_ID", "ses_index")], by = "Poly_ID")
sofia_income_noise$ses_index <- -sofia_income_noise$ses_index

# Plot data
qtm(sofia_income_noise, "noise_popw")
qtm(sofia_income_noise, "ses_index")

# Spatial weights matrix
centroids <- st_centroid(st_geometry(sofia_income_noise))
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
isolated_areas <- which(card(wm_d) == 0) # Identify areas with zero neighbors
sofia_income_noise <- sofia_income_noise[-isolated_areas, ] # Remove isolated areas from the data
centroids <- centroids[-isolated_areas, ]
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
# Row standardized weights matrix 
rswm_d <- nb2listw(wm_d, style = "W", zero.policy = TRUE)
summary(rswm_d)

### Analysis for Noise
# Bivariate Moran I scatter plot
sofia_income_noise$std_ses <- standardize(sofia_income_noise$ses_index)
sofia_income_noise$std_noise <- standardize(sofia_income_noise$noise_popw)
# Lag variable: 
sofia_income_noise$wsr_noise <- lag.listw(rswm_d,sofia_income_noise$std_noise)

# Plot: 
theme_set(theme_bw())
p3 <- ggplot(data = sofia_income_noise, aes(x=std_ses,y=wsr_noise)) +
  geom_point(color = "grey", alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty =2) +
  ggtitle("SES and lagged noise") + ylab("Standard lagged noise") + xlab("Standard SES index")

### Global Moran I (bivariate)
bivariate_moran <- moran_bv(sofia_income_noise$std_ses, sofia_income_noise$std_noise, rswm_d, nsim = 500)
print(bivariate_moran)
bivariate_moran$t0
summary(lm(sofia_income_noise$wsr_noise ~ sofia_income_noise$std_ses))

### Local Moran I (bivariate)
localMI <- localmoran_bv(sofia_income_noise$std_ses, sofia_income_noise$std_noise, rswm_d, nsim = 100)
head(localMI)

# Plot results
sofia_income_noise.localMI <- cbind(sofia_income_noise, localMI)

# Moran I
localMI.map <- tm_shape(sofia_income_noise.localMI) +
  tm_fill(col = "Ibvi",
          style = "pretty",
          title = "local moran statistics")+
  tm_borders(alpha = 0.5)
# p values
pvalue.map <- tm_shape(sofia_income_noise.localMI) +
  tm_fill(col = "Pr.z....E.Ibvi..",
          breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
          palette = "-Blues",
          title = "local Moran's I p-values")+
  tm_borders(alpha = 0.5)

tmap_arrange(localMI.map, pvalue.map, asp = 1, ncol = 2)

### LISA cluster map
# Identifying the LISA clusters
quadrant <- vector(mode = "numeric", length = nrow(localMI))
signif <- 0.05 # significance level

# Assign quadrant values (HH, LL, HL, LH)
quadrant[sofia_income_noise.localMI$std_ses > 0 & sofia_income_noise.localMI$wsr_noise > 0] <- "High SES-High Noise"      
quadrant[sofia_income_noise.localMI$std_ses < 0 & sofia_income_noise.localMI$wsr_noise < 0] <- "Low SES-Low Noise"      
quadrant[sofia_income_noise.localMI$std_ses < 0 & sofia_income_noise.localMI$wsr_noise > 0] <- "Low SES-High Noise"
quadrant[sofia_income_noise.localMI$std_ses > 0 & sofia_income_noise.localMI$wsr_noise < 0] <- "High SES-Low Noise"
quadrant[localMI[,5]>signif] <- "NS"
sofia_income_noise.localMI$quadrant <- quadrant 

# Merge back to full polygon
sofia_income_noise.localMI$geometry <- NULL
sofia_income_noise <- merge(sofia_noise, sofia_income_noise.localMI[c("Poly_ID", "quadrant")], by = "Poly_ID", all.x = T)
sofia_income_noise <- sofia_income_noise %>% mutate(quadrant = replace_na(quadrant, "NS"))
sofia_income_noise$quadrant <- factor(sofia_income_noise$quadrant, levels = c("Low SES-High Noise", "High SES-High Noise", "High SES-Low Noise", "Low SES-Low Noise", "NS"))

# Make the plot
theme_set(theme_bw())
w5 <- ggplot(sofia_income_noise) + geom_sf(aes(fill = quadrant), size = 0.1) + 
  scale_fill_manual(name = "", values = c("High SES-High Noise" = "#fdae61",  
                                          "Low SES-Low Noise" = "#abd9e9", 
                                          "Low SES-High Noise" = "#d7191c", 
                                          "High SES-Low Noise" = "#2c7bb6", 
                                          "NS" = "#ffffff"), guide = guide_legend(nrow = 2, keywidth = 1.5, keyheight = 1.5)) +
  ggtitle("Road-traffic noise exposure") + 
  theme(legend.position = "bottom", legend.text = element_text(size = 12), 
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        panel.grid = element_blank()) + 
  annotation_scale(location = "bl", width_hint = 0.2) +  # Add scale bar at bottom left
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering())

ggsave("Baseline_BOD/Sofia_spatial_analysis_Noise.png")


#-----------------------------------
### Green space levels
sofia_green <- st_read("Baseline_BOD/Clean_data/sofia_green_poly.geojson")
sofia_income_green <- merge(sofia_green, sofia_ses_poly[c("Poly_ID", "ses_index")], by = "Poly_ID")
sofia_income_green$ses_index <- -sofia_income_green$ses_index

# Plot data
qtm(sofia_income_green, "ndvi")
qtm(sofia_income_green, "ses_index")

# Spatial weights matrix
centroids <- st_centroid(st_geometry(sofia_income_green))
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
isolated_areas <- which(card(wm_d) == 0) # Identify areas with zero neighbors
sofia_income_green <- sofia_income_green[-isolated_areas, ] # Remove isolated areas from the data
centroids <- centroids[-isolated_areas, ]
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
# Row standardized weights matrix 
rswm_d <- nb2listw(wm_d, style = "W", zero.policy = TRUE)
summary(rswm_d)

### Analysis for green space
# Bivariate Moran I scatter plot
sofia_income_green$std_ses <- standardize(sofia_income_green$ses_index)
sofia_income_green$std_ndvi <- standardize(sofia_income_green$ndvi)
# Lag variable: 
sofia_income_green$wsr_ndvi <- lag.listw(rswm_d, sofia_income_green$std_ndvi)

# Plot: 
theme_set(theme_bw())
p4 <- ggplot(data = sofia_income_green, aes(x=std_ses, y=wsr_ndvi)) +
  geom_point(color = "grey", alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty =2) +
  ggtitle("SES and lagged NDVI") + ylab("Standard lagged NDVI") + xlab("Standard SES index")

### Global Moran I (bivariate)
bivariate_moran <- moran_bv(sofia_income_green$std_ses, sofia_income_green$std_ndvi, rswm_d, nsim = 500)
print(bivariate_moran)
bivariate_moran$t0
summary(lm(sofia_income_green$wsr_ndvi ~ sofia_income_green$std_ses))

### Local Moran I (bivariate)
localMI <- localmoran_bv(sofia_income_green$std_ses, sofia_income_green$std_ndvi, rswm_d, nsim = 100)
head(localMI)

# Plot results
sofia_income_green.localMI <- cbind(sofia_income_green, localMI)

# Moran I
localMI.map <- tm_shape(sofia_income_green.localMI) +
  tm_fill(col = "Ibvi",
          style = "pretty",
          title = "local moran statistics")+
  tm_borders(alpha = 0.5)
# p values
pvalue.map <- tm_shape(sofia_income_green.localMI) +
  tm_fill(col = "Pr.z....E.Ibvi..",
          breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
          palette = "-Blues",
          title = "local Moran's I p-values")+
  tm_borders(alpha = 0.5)

tmap_arrange(localMI.map, pvalue.map, asp = 1, ncol = 2)

### LISA cluster map
# Identifying the LISA clusters
quadrant <- vector(mode = "numeric", length = nrow(localMI))
signif <- 0.05 # significance level

# Assign quadrant values (HH, LL, HL, LH)
quadrant[sofia_income_green.localMI$std_ses > 0 & sofia_income_green.localMI$wsr_ndvi > 0] <- "High SES-High NDVI"      
quadrant[sofia_income_green.localMI$std_ses < 0 & sofia_income_green.localMI$wsr_ndvi < 0] <- "Low SES-Low NDVI"      
quadrant[sofia_income_green.localMI$std_ses < 0 & sofia_income_green.localMI$wsr_ndvi > 0] <- "Low SES-High NDVI"
quadrant[sofia_income_green.localMI$std_ses > 0 & sofia_income_green.localMI$wsr_ndvi < 0] <- "High SES-Low NDVI"
quadrant[localMI[,5]>signif] <- "NS"
sofia_income_green.localMI$quadrant <- quadrant 

# Merge back to full polygon
sofia_income_green.localMI$geometry <- NULL
sofia_income_green <- merge(sofia_green, sofia_income_green.localMI[c("Poly_ID", "quadrant")], by = "Poly_ID", all.x = T)
sofia_income_green <- sofia_income_green %>% mutate(quadrant = replace_na(quadrant, "NS"))
sofia_income_green$quadrant <- factor(sofia_income_green$quadrant, levels = c("Low SES-Low NDVI", "High SES-Low NDVI", "High SES-High NDVI", "Low SES-High NDVI", "NS"))

# Make the plot
theme_set(theme_bw())
w7 <- ggplot(sofia_income_green) + geom_sf(aes(fill = quadrant), size = 0.1) + 
  scale_fill_manual(name = "", values = c("High SES-Low NDVI" = "#fdae61",  
                                          "Low SES-High NDVI" = "#abd9e9", 
                                          "Low SES-Low NDVI" = "#d7191c", 
                                          "High SES-High NDVI" = "#2c7bb6", 
                                          "NS" = "#ffffff"), guide = guide_legend(nrow = 2, keywidth = 1.5, keyheight = 1.5)) +
  ggtitle("NDVI exposure") + 
  theme(legend.position = "bottom", legend.text = element_text(size = 12), 
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        panel.grid = element_blank()) + 
  annotation_scale(location = "bl", width_hint = 0.2) +  # Add scale bar at bottom left
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering())

ggsave("Baseline_BOD/Sofia_spatial_analysis_NDVI.png")


#-----------------------------------
### Green space levels
sofia_green <- st_read("Baseline_BOD/Clean_data/sofia_green_poly.geojson")
sofia_income_green <- merge(sofia_green, sofia_ses_poly[c("Poly_ID", "ses_index")], by = "Poly_ID")
sofia_income_green$ses_index <- -sofia_income_green$ses_index

# Plot data
qtm(sofia_income_green, "per_green")
qtm(sofia_income_green, "ses_index")

# Spatial weights matrix
centroids <- st_centroid(st_geometry(sofia_income_green))
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
isolated_areas <- which(card(wm_d) == 0) # Identify areas with zero neighbors
sofia_income_green <- sofia_income_green[-isolated_areas, ] # Remove isolated areas from the data
centroids <- centroids[-isolated_areas, ]
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
# Row standardized weights matrix 
rswm_d <- nb2listw(wm_d, style = "W", zero.policy = TRUE)
summary(rswm_d)

### Analysis for green space
# Bivariate Moran I scatter plot
sofia_income_green$std_ses <- standardize(sofia_income_green$ses_index)
sofia_income_green$std_green <- standardize(sofia_income_green$per_green)
# Lag variable: 
sofia_income_green$wsr_green <- lag.listw(rswm_d, sofia_income_green$std_green)

# Plot: 
theme_set(theme_bw())
p5 <- ggplot(data = sofia_income_green, aes(x=std_ses, y=wsr_green)) +
  geom_point(color = "grey", alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty =2) +
  ggtitle("SES and lagged %GA") + ylab("Standard lagged %GA") + xlab("Standard SES index")

### Global Moran I (bivariate)
bivariate_moran <- moran_bv(sofia_income_green$std_ses, sofia_income_green$std_green, rswm_d, nsim = 500)
print(bivariate_moran)
bivariate_moran$t0
summary(lm(sofia_income_green$wsr_green ~ sofia_income_green$std_ses))

### Local Moran I (bivariate)
localMI <- localmoran_bv(sofia_income_green$std_ses, sofia_income_green$std_green, rswm_d, nsim = 100)
head(localMI)

# Plot results
sofia_income_green.localMI <- cbind(sofia_income_green, localMI)

# Moran I
localMI.map <- tm_shape(sofia_income_green.localMI) +
  tm_fill(col = "Ibvi",
          style = "pretty",
          title = "local moran statistics")+
  tm_borders(alpha = 0.5)
# p values
pvalue.map <- tm_shape(sofia_income_green.localMI) +
  tm_fill(col = "Pr.z....E.Ibvi..",
          breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
          palette = "-Blues",
          title = "local Moran's I p-values")+
  tm_borders(alpha = 0.5)

tmap_arrange(localMI.map, pvalue.map, asp = 1, ncol = 2)

### LISA cluster map
# Identifying the LISA clusters
quadrant <- vector(mode = "numeric", length = nrow(localMI))
signif <- 0.05 # significance level

# Assign quadrant values (HH, LL, HL, LH)
quadrant[sofia_income_green.localMI$std_ses > 0 & sofia_income_green.localMI$wsr_green > 0] <- "High SES-High %GA"      
quadrant[sofia_income_green.localMI$std_ses < 0 & sofia_income_green.localMI$wsr_green < 0] <- "Low SES-Low %GA"      
quadrant[sofia_income_green.localMI$std_ses < 0 & sofia_income_green.localMI$wsr_green > 0] <- "Low SES-High %GA"
quadrant[sofia_income_green.localMI$std_ses > 0 & sofia_income_green.localMI$wsr_green < 0] <- "High SES-Low %GA"
quadrant[localMI[,5]>signif] <- "NS"
sofia_income_green.localMI$quadrant <- quadrant 

# Merge back to full polygon
sofia_income_green.localMI$geometry <- NULL
sofia_income_green <- merge(sofia_green, sofia_income_green.localMI[c("Poly_ID", "quadrant")], by = "Poly_ID", all.x = T)
sofia_income_green <- sofia_income_green %>% mutate(quadrant = replace_na(quadrant, "NS"))
sofia_income_green$quadrant <- factor(sofia_income_green$quadrant, levels = c("Low SES-Low %GA", "High SES-Low %GA", "High SES-High %GA", "Low SES-High %GA", "NS"))

# Make the plot
theme_set(theme_bw())
w9 <- ggplot(sofia_income_green) + geom_sf(aes(fill = quadrant), size = 0.1) + 
  scale_fill_manual(name = "", values = c("High SES-Low %GA" = "#fdae61",  
                                          "Low SES-High %GA" = "#abd9e9", 
                                          "Low SES-Low %GA" = "#d7191c", 
                                          "High SES-High %GA" = "#2c7bb6", 
                                          "NS" = "#ffffff"), guide = guide_legend(nrow = 2, keywidth = 1.5, keyheight = 1.5)) +
  ggtitle("%GA exposure") + 
  theme(legend.position = "bottom", legend.text = element_text(size = 12), 
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        panel.grid = element_blank()) + 
  annotation_scale(location = "bl", width_hint = 0.2) +  # Add scale bar at bottom left
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering())

ggsave("Baseline_BOD/Sofia_spatial_analysis_%GA.png")


#-----------------------------------
### UHI effect
sofia_uhi <- st_read("Baseline_BOD/Clean_data/sofia_uhi_poly.geojson")
sofia_income_uhi <- merge(sofia_uhi, sofia_ses_poly[c("Poly_ID", "ses_index")], by = "Poly_ID")
sofia_income_uhi$ses_index <- -sofia_income_uhi$ses_index

# Plot data
qtm(sofia_income_uhi, "uhi_summer")
qtm(sofia_income_uhi, "ses_index")

# Spatial weights matrix
centroids <- st_centroid(st_geometry(sofia_income_uhi))
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
isolated_areas <- which(card(wm_d) == 0) # Identify areas with zero neighbors
sofia_income_uhi <- sofia_income_uhi[-isolated_areas, ] # Remove isolated areas from the data
centroids <- centroids[-isolated_areas, ]
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
# Row standardized weights matrix 
rswm_d <- nb2listw(wm_d, style = "W", zero.policy = TRUE)
summary(rswm_d)

### Analysis for green space
# Bivariate Moran I scatter plot
sofia_income_uhi$std_ses <- standardize(sofia_income_uhi$ses_index)
sofia_income_uhi$std_uhi <- standardize(sofia_income_uhi$uhi_summer)
# Lag variable: 
sofia_income_uhi$wsr_uhi <- lag.listw(rswm_d, sofia_income_uhi$std_uhi)

# Plot: 
theme_set(theme_bw())
p6 <- ggplot(data = sofia_income_uhi, aes(x=std_ses, y=wsr_uhi)) +
  geom_point(color = "grey", alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty =2) +
  ggtitle("SES and lagged UHI") + ylab("Standard lagged UHI") + xlab("Standard SES index")

### Global Moran I (bivariate)
bivariate_moran <- moran_bv(sofia_income_uhi$std_ses, sofia_income_uhi$std_uhi, rswm_d, nsim = 500)
print(bivariate_moran)
bivariate_moran$t0
summary(lm(sofia_income_uhi$wsr_uhi ~ sofia_income_uhi$std_ses))

### Local Moran I (bivariate)
localMI <- localmoran_bv(sofia_income_uhi$std_ses, sofia_income_uhi$std_uhi, rswm_d, nsim = 100)
head(localMI)

# Plot results
sofia_income_uhi.localMI <- cbind(sofia_income_uhi, localMI)

# Moran I
localMI.map <- tm_shape(sofia_income_uhi.localMI) +
  tm_fill(col = "Ibvi",
          style = "pretty",
          title = "local moran statistics")+
  tm_borders(alpha = 0.5)
# p values
pvalue.map <- tm_shape(sofia_income_uhi.localMI) +
  tm_fill(col = "Pr.z....E.Ibvi..",
          breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
          palette = "-Blues",
          title = "local Moran's I p-values")+
  tm_borders(alpha = 0.5)

tmap_arrange(localMI.map, pvalue.map, asp = 1, ncol = 2)

### LISA cluster map
# Identifying the LISA clusters
quadrant <- vector(mode = "numeric", length = nrow(localMI))
signif <- 0.05 # significance level

# Assign quadrant values (HH, LL, HL, LH)
quadrant[sofia_income_uhi.localMI$std_ses > 0 & sofia_income_uhi.localMI$wsr_uhi > 0] <- "High SES-High UHI"      
quadrant[sofia_income_uhi.localMI$std_ses < 0 & sofia_income_uhi.localMI$wsr_uhi < 0] <- "Low SES-Low UHI"      
quadrant[sofia_income_uhi.localMI$std_ses < 0 & sofia_income_uhi.localMI$wsr_uhi > 0] <- "Low SES-High UHI"
quadrant[sofia_income_uhi.localMI$std_ses > 0 & sofia_income_uhi.localMI$wsr_uhi < 0] <- "High SES-Low UHI"
quadrant[localMI[,5]>signif] <- "NS"
sofia_income_uhi.localMI$quadrant <- quadrant 

# Merge back to full polygon
sofia_income_uhi.localMI$geometry <- NULL
sofia_income_uhi <- merge(sofia_uhi, sofia_income_uhi.localMI[c("Poly_ID", "quadrant")], by = "Poly_ID", all.x = T)
sofia_income_uhi <- sofia_income_uhi %>% mutate(quadrant = replace_na(quadrant, "NS"))
sofia_income_uhi$quadrant <- factor(sofia_income_uhi$quadrant, levels = c("Low SES-High UHI", "High SES-High UHI", "High SES-Low UHI", "Low SES-Low UHI", "NS"))

# Make the plot
theme_set(theme_bw())
w11 <- ggplot(sofia_income_uhi) + geom_sf(aes(fill = quadrant), size = 0.1) + 
  scale_fill_manual(name = "", values = c("High SES-High UHI" = "#fdae61",  
                                          "Low SES-Low UHI" = "#abd9e9", 
                                          "Low SES-High UHI" = "#d7191c", 
                                          "High SES-Low UHI" = "#2c7bb6", 
                                          "NS" = "#ffffff"), guide = guide_legend(nrow = 2, keywidth = 1.5, keyheight = 1.5)) +
  ggtitle("UHI exposure") + 
  theme(legend.position = "bottom", legend.text = element_text(size = 12), 
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        panel.grid = element_blank()) + 
  annotation_scale(location = "bl", width_hint = 0.2) +  # Add scale bar at bottom left
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering())

ggsave("Baseline_BOD/Sofia_spatial_analysis_UHI.png")

### Plot with all the exposures
library(ggpubr)
ggarrange(w1, w3, w5, w7, w9, w11,  ncol = 3, nrow = 2)
ggsave("Baseline_BOD/Sofia_spatial_analysis_polygons_exposures.png")

# Bivariate Moran I biplot
ggarrange(p1, p2, p3, p4, p5, p6,  ncol = 2, nrow = 3)
ggsave("Baseline_BOD/Sofia_bivariateMoran_biplot_exposures.png")


#-----------------------------------
### PM2.5 mortality
sofia_airp <- st_read("Baseline_BOD/Results/res_pm25_nat_mort.geojson")
sofia_airp <- st_transform(sofia_airp, crs = 3035)
sofia_income_airp <- merge(sofia_airp, sofia_ses_poly[c("Poly_ID", "ses_index")], by = "Poly_ID")
sofia_income_airp$ses_index <- -sofia_income_airp$ses_index

# Plot data
qtm(sofia_income_airp, "att_rate")
qtm(sofia_income_airp, "ses_index")

# Spatial weights matrix
centroids <- st_centroid(st_geometry(sofia_income_airp))
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
isolated_areas <- which(card(wm_d) == 0) # Identify areas with zero neighbors
sofia_income_airp <- sofia_income_airp[-isolated_areas, ] # Remove isolated areas from the data
centroids <- centroids[-isolated_areas, ]
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
# Row standardized weights matrix 
rswm_d <- nb2listw(wm_d, style = "W", zero.policy = TRUE)
summary(rswm_d)

### Analysis for PM2.5
# Bivariate Moran I scatter plot
sofia_income_airp$std_ses <- standardize(sofia_income_airp$ses_index)
sofia_income_airp$std_rate <- standardize(sofia_income_airp$att_rate)
# Lag variable: 
sofia_income_airp$wsr_rate <- lag.listw(rswm_d,sofia_income_airp$std_rate)

# Plot: 
theme_set(theme_bw())
p1 <- ggplot(data = sofia_income_airp, aes(x=std_ses,y=wsr_rate)) +
  geom_point(color = "grey", alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty =2) +
  ggtitle(expression("SES and lagged PM"[2.5]*" mortality")) + ylab(expression("Standard lagged PM"[2.5]*" mortality")) + xlab("Standard SES index")

### Global Moran I (bivariate)
bivariate_moran <- moran_bv(sofia_income_airp$std_ses, sofia_income_airp$std_rate, rswm_d, nsim = 500)
print(bivariate_moran)
bivariate_moran$t0
summary(lm(sofia_income_airp$wsr_rate ~ sofia_income_airp$std_ses))

### Local Moran I (bivariate)
localMI <- localmoran_bv(sofia_income_airp$std_ses, sofia_income_airp$std_rate, rswm_d, nsim = 100)
head(localMI)

# Plot results
sofia_income_airp.localMI <- cbind(sofia_income_airp, localMI)

# Moran I
localMI.map <- tm_shape(sofia_income_airp.localMI) +
  tm_fill(col = "Ibvi",
          style = "pretty",
          title = "local moran statistics")+
  tm_borders(alpha = 0.5)
# p values
pvalue.map <- tm_shape(sofia_income_airp.localMI) +
  tm_fill(col = "Pr.z....E.Ibvi..",
          breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
          palette = "-Blues",
          title = "local Moran's I p-values")+
  tm_borders(alpha = 0.5)

tmap_arrange(localMI.map, pvalue.map, asp = 1, ncol = 2)

### LISA cluster map
# Identifying the LISA clusters
quadrant <- vector(mode = "numeric", length = nrow(localMI))
signif <- 0.05 # significance level

# Assign quadrant values (HH, LL, HL, LH)
quadrant[sofia_income_airp.localMI$std_ses > 0 & sofia_income_airp.localMI$wsr_rate > 0] <- "High SES-High Mortality"      
quadrant[sofia_income_airp.localMI$std_ses < 0 & sofia_income_airp.localMI$wsr_rate < 0] <- "Low SES-Low Mortality"      
quadrant[sofia_income_airp.localMI$std_ses < 0 & sofia_income_airp.localMI$wsr_rate > 0] <- "Low SES-High Mortality"
quadrant[sofia_income_airp.localMI$std_ses > 0 & sofia_income_airp.localMI$wsr_rate < 0] <- "High SES-Low Mortality"
quadrant[localMI[,5]>signif] <- "NS"
sofia_income_airp.localMI$quadrant <- quadrant 

# Merge back to full polygon
sofia_income_airp.localMI$geometry <- NULL
sofia_income_airp <- merge(sofia_airp, sofia_income_airp.localMI[c("Poly_ID", "quadrant")], by = "Poly_ID", all.x = T)
sofia_income_airp <- sofia_income_airp %>% mutate(quadrant = replace_na(quadrant, "NS"))
sofia_income_airp$quadrant <- factor(sofia_income_airp$quadrant, levels = c("Low SES-High Mortality", "High SES-High Mortality", "High SES-Low Mortality", "Low SES-Low Mortality", "NS"))

# Make the plot
theme_set(theme_bw())
w2 <- ggplot(sofia_income_airp) + geom_sf(aes(fill = quadrant), size = 0.1) + 
  scale_fill_manual(name = "", values = c("High SES-High Mortality" = "#fdae61",  
                                          "Low SES-Low Mortality" = "#abd9e9", 
                                          "Low SES-High Mortality" = "#d7191c", 
                                          "High SES-Low Mortality" = "#2c7bb6", 
                                          "NS" = "#ffffff"), guide = guide_legend(nrow = 2, keywidth = 1.5, keyheight = 1.5)) +
  ggtitle(expression("PM"[2.5]*" mortality")) + 
  theme(legend.position = "bottom", legend.text = element_text(size = 11), 
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        panel.grid = element_blank()) + 
  annotation_scale(location = "bl", width_hint = 0.2) +  # Add scale bar at bottom left
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering())

ggsave("Baseline_BOD/Sofia_spatial_analysis_PM2.5_mortality.png")


#-----------------------------------
### NO2 mortality
sofia_airp <- st_read("Baseline_BOD/Results/res_no2_nat_mort.geojson")
sofia_airp <- st_transform(sofia_airp, crs = 3035)
sofia_income_airp <- merge(sofia_airp, sofia_ses_poly[c("Poly_ID", "ses_index")], by = "Poly_ID")
sofia_income_airp$ses_index <- -sofia_income_airp$ses_index

# Plot data
qtm(sofia_income_airp, "att_rate")
qtm(sofia_income_airp, "ses_index")

# Spatial weights matrix
centroids <- st_centroid(st_geometry(sofia_income_airp))
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
isolated_areas <- which(card(wm_d) == 0) # Identify areas with zero neighbors
sofia_income_airp <- sofia_income_airp[-isolated_areas, ] # Remove isolated areas from the data
centroids <- centroids[-isolated_areas, ]
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
# Row standardized weights matrix 
rswm_d <- nb2listw(wm_d, style = "W", zero.policy = TRUE)
summary(rswm_d)

### Analysis for NO2
# Bivariate Moran I scatter plot
sofia_income_airp$std_ses <- standardize(sofia_income_airp$ses_index)
sofia_income_airp$std_rate <- standardize(sofia_income_airp$att_rate)
# Lag variable: 
sofia_income_airp$wsr_rate <- lag.listw(rswm_d,sofia_income_airp$std_rate)

# Plot: 
theme_set(theme_bw())
p2 <- ggplot(data = sofia_income_airp, aes(x=std_ses,y=wsr_rate)) +
  geom_point(color = "grey", alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty =2) +
  ggtitle(expression("SES and lagged NO"[2]*" mortality")) + ylab(expression("Standard lagged NO"[2]*" mortality")) + xlab("Standard SES index")

### Global Moran I (bivariate)
bivariate_moran <- moran_bv(sofia_income_airp$std_ses, sofia_income_airp$std_rate, rswm_d, nsim = 500)
print(bivariate_moran)
bivariate_moran$t0
summary(lm(sofia_income_airp$wsr_rate ~ sofia_income_airp$std_ses))

### Local Moran I (bivariate)
localMI <- localmoran_bv(sofia_income_airp$std_ses, sofia_income_airp$std_rate, rswm_d, nsim = 100)
head(localMI)

# Plot results
sofia_income_airp.localMI <- cbind(sofia_income_airp, localMI)

# Moran I
localMI.map <- tm_shape(sofia_income_airp.localMI) +
  tm_fill(col = "Ibvi",
          style = "pretty",
          title = "local moran statistics")+
  tm_borders(alpha = 0.5)
# p values
pvalue.map <- tm_shape(sofia_income_airp.localMI) +
  tm_fill(col = "Pr.z....E.Ibvi..",
          breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
          palette = "-Blues",
          title = "local Moran's I p-values")+
  tm_borders(alpha = 0.5)

tmap_arrange(localMI.map, pvalue.map, asp = 1, ncol = 2)

### LISA cluster map
# Identifying the LISA clusters
quadrant <- vector(mode = "numeric", length = nrow(localMI))
signif <- 0.05 # significance level

# Assign quadrant values (HH, LL, HL, LH)
quadrant[sofia_income_airp.localMI$std_ses > 0 & sofia_income_airp.localMI$wsr_rate > 0] <- "High SES-High Mortality"      
quadrant[sofia_income_airp.localMI$std_ses < 0 & sofia_income_airp.localMI$wsr_rate < 0] <- "Low SES-Low Mortality"      
quadrant[sofia_income_airp.localMI$std_ses < 0 & sofia_income_airp.localMI$wsr_rate > 0] <- "Low SES-High Mortality"
quadrant[sofia_income_airp.localMI$std_ses > 0 & sofia_income_airp.localMI$wsr_rate < 0] <- "High SES-Low Mortality"
quadrant[localMI[,5]>signif] <- "NS"
sofia_income_airp.localMI$quadrant <- quadrant 

# Merge back to full polygon
sofia_income_airp.localMI$geometry <- NULL
sofia_income_airp <- merge(sofia_airp, sofia_income_airp.localMI[c("Poly_ID", "quadrant")], by = "Poly_ID", all.x = T)
sofia_income_airp <- sofia_income_airp %>% mutate(quadrant = replace_na(quadrant, "NS"))
sofia_income_airp$quadrant <- factor(sofia_income_airp$quadrant, levels = c("Low SES-High Mortality", "High SES-High Mortality", "High SES-Low Mortality", "Low SES-Low Mortality", "NS"))

# Make the plot
theme_set(theme_bw())
w4 <- ggplot(sofia_income_airp) + geom_sf(aes(fill = quadrant), size = 0.1) + 
  scale_fill_manual(name = "", values = c("High SES-High Mortality" = "#fdae61",  
                                          "Low SES-Low Mortality" = "#abd9e9", 
                                          "Low SES-High Mortality" = "#d7191c", 
                                          "High SES-Low Mortality" = "#2c7bb6", 
                                          "NS" = "#ffffff"), guide = guide_legend(nrow = 2, keywidth = 1.5, keyheight = 1.5)) +
  ggtitle(expression("NO"[2]*" mortality")) + 
  theme(legend.position = "bottom", legend.text = element_text(size = 11), 
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        panel.grid = element_blank()) + 
  annotation_scale(location = "bl", width_hint = 0.2) +  # Add scale bar at bottom left
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering())

ggsave("Baseline_BOD/Sofia_spatial_analysis_NO2_mortality.png")

#-----------------------------------
### PM2.5 and NO2 mortality
sofia_airp <- st_read("Baseline_BOD/Results/res_pm25_nat_mort_multiexp.geojson")
sofia_airp <- st_transform(sofia_airp, crs = 3035)
sofia_income_airp <- merge(sofia_airp, sofia_ses_poly[c("Poly_ID", "ses_index")], by = "Poly_ID")
sofia_income_airp$ses_index <- -sofia_income_airp$ses_index

# Plot data
qtm(sofia_income_airp, "att_rate")
qtm(sofia_income_airp, "ses_index")

# Spatial weights matrix
centroids <- st_centroid(st_geometry(sofia_income_airp))
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
isolated_areas <- which(card(wm_d) == 0) # Identify areas with zero neighbors
sofia_income_airp <- sofia_income_airp[-isolated_areas, ] # Remove isolated areas from the data
centroids <- centroids[-isolated_areas, ]
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
# Row standardized weights matrix 
rswm_d <- nb2listw(wm_d, style = "W", zero.policy = TRUE)
summary(rswm_d)

### Analysis for NO2
# Bivariate Moran I scatter plot
sofia_income_airp$std_ses <- standardize(sofia_income_airp$ses_index)
sofia_income_airp$std_rate <- standardize(sofia_income_airp$att_rate)
# Lag variable: 
sofia_income_airp$wsr_rate <- lag.listw(rswm_d,sofia_income_airp$std_rate)

# Plot: 
theme_set(theme_bw())
p3 <- ggplot(data = sofia_income_airp, aes(x=std_ses,y=wsr_rate)) +
  geom_point(color = "grey", alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty =2) +
  ggtitle(expression("SES and lagged PM"[2.5]*" and NO"[2]*" mortality")) + ylab(expression("Standard lagged PM"[2.5]*" and NO"[2]*" mortality")) + xlab("Standard SES index")

### Global Moran I (bivariate)
bivariate_moran <- moran_bv(sofia_income_airp$std_ses, sofia_income_airp$std_rate, rswm_d, nsim = 500)
print(bivariate_moran)
bivariate_moran$t0
summary(lm(sofia_income_airp$wsr_rate ~ sofia_income_airp$std_ses))

### Local Moran I (bivariate)
localMI <- localmoran_bv(sofia_income_airp$std_ses, sofia_income_airp$std_rate, rswm_d, nsim = 100)
head(localMI)

# Plot results
sofia_income_airp.localMI <- cbind(sofia_income_airp, localMI)

# Moran I
localMI.map <- tm_shape(sofia_income_airp.localMI) +
  tm_fill(col = "Ibvi",
          style = "pretty",
          title = "local moran statistics")+
  tm_borders(alpha = 0.5)
# p values
pvalue.map <- tm_shape(sofia_income_airp.localMI) +
  tm_fill(col = "Pr.z....E.Ibvi..",
          breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
          palette = "-Blues",
          title = "local Moran's I p-values")+
  tm_borders(alpha = 0.5)

tmap_arrange(localMI.map, pvalue.map, asp = 1, ncol = 2)

### LISA cluster map
# Identifying the LISA clusters
quadrant <- vector(mode = "numeric", length = nrow(localMI))
signif <- 0.05 # significance level

# Assign quadrant values (HH, LL, HL, LH)
quadrant[sofia_income_airp.localMI$std_ses > 0 & sofia_income_airp.localMI$wsr_rate > 0] <- "High SES-High Mortality"      
quadrant[sofia_income_airp.localMI$std_ses < 0 & sofia_income_airp.localMI$wsr_rate < 0] <- "Low SES-Low Mortality"      
quadrant[sofia_income_airp.localMI$std_ses < 0 & sofia_income_airp.localMI$wsr_rate > 0] <- "Low SES-High Mortality"
quadrant[sofia_income_airp.localMI$std_ses > 0 & sofia_income_airp.localMI$wsr_rate < 0] <- "High SES-Low Mortality"
quadrant[localMI[,5]>signif] <- "NS"
sofia_income_airp.localMI$quadrant <- quadrant 

# Merge back to full polygon
sofia_income_airp.localMI$geometry <- NULL
sofia_income_airp <- merge(sofia_airp, sofia_income_airp.localMI[c("Poly_ID", "quadrant")], by = "Poly_ID", all.x = T)
sofia_income_airp <- sofia_income_airp %>% mutate(quadrant = replace_na(quadrant, "NS"))
sofia_income_airp$quadrant <- factor(sofia_income_airp$quadrant, levels = c("Low SES-High Mortality", "High SES-High Mortality", "High SES-Low Mortality", "Low SES-Low Mortality", "NS"))

# Make the plot
theme_set(theme_bw())
w6 <- ggplot(sofia_income_airp) + geom_sf(aes(fill = quadrant), size = 0.1) + 
  scale_fill_manual(name = "", values = c("High SES-High Mortality" = "#fdae61",  
                                          "Low SES-Low Mortality" = "#abd9e9", 
                                          "Low SES-High Mortality" = "#d7191c", 
                                          "High SES-Low Mortality" = "#2c7bb6", 
                                          "NS" = "#ffffff"), guide = guide_legend(nrow = 2, keywidth = 1.5, keyheight = 1.5)) +
  ggtitle(expression("PM"[2.5]*" and NO"[2]*" mortality")) + 
  theme(legend.position = "bottom", legend.text = element_text(size = 11), 
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        panel.grid = element_blank()) + 
  annotation_scale(location = "bl", width_hint = 0.2) +  # Add scale bar at bottom left
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering())

ggsave("Baseline_BOD/Sofia_spatial_analysis_PM25_NO2_mortality.png")


#-----------------------------------
### Noise mortality
sofia_noise <- st_read("Baseline_BOD/Results/res_noise_nat_mort.geojson")
sofia_noise <- st_transform(sofia_noise, crs = 3035)
sofia_income_noise <- merge(sofia_noise, sofia_ses_poly[c("Poly_ID", "ses_index")], by = "Poly_ID")
sofia_income_noise$ses_index <- -sofia_income_noise$ses_index

# Plot data
qtm(sofia_income_noise, "att_rate")
qtm(sofia_income_noise, "ses_index")

# Spatial weights matrix
centroids <- st_centroid(st_geometry(sofia_income_noise))
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
isolated_areas <- which(card(wm_d) == 0) # Identify areas with zero neighbors
sofia_income_noise <- sofia_income_noise[-isolated_areas, ] # Remove isolated areas from the data
centroids <- centroids[-isolated_areas, ]
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
# Row standardized weights matrix 
rswm_d <- nb2listw(wm_d, style = "W", zero.policy = TRUE)
summary(rswm_d)

### Analysis for noise
# Bivariate Moran I scatter plot
sofia_income_noise$std_ses <- standardize(sofia_income_noise$ses_index)
sofia_income_noise$std_rate <- standardize(sofia_income_noise$att_rate)
# Lag variable: 
sofia_income_noise$wsr_rate <- lag.listw(rswm_d,sofia_income_noise$std_rate)

# Plot: 
theme_set(theme_bw())
p4 <- ggplot(data = sofia_income_noise, aes(x=std_ses,y=wsr_rate)) +
  geom_point(color = "grey", alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty =2) +
  ggtitle("SES and lagged noise mortality") + ylab("Standard lagged noise mortality") + xlab("Standard SES index")

### Global Moran I (bivariate)
bivariate_moran <- moran_bv(sofia_income_noise$std_ses, sofia_income_noise$std_rate, rswm_d, nsim = 500)
print(bivariate_moran)
bivariate_moran$t0
summary(lm(sofia_income_noise$wsr_rate ~ sofia_income_noise$std_ses))

### Local Moran I (bivariate)
localMI <- localmoran_bv(sofia_income_noise$std_ses, sofia_income_noise$std_rate, rswm_d, nsim = 100)
head(localMI)

# Plot results
sofia_income_noise.localMI <- cbind(sofia_income_noise, localMI)

# Moran I
localMI.map <- tm_shape(sofia_income_noise.localMI) +
  tm_fill(col = "Ibvi",
          style = "pretty",
          title = "local moran statistics")+
  tm_borders(alpha = 0.5)
# p values
pvalue.map <- tm_shape(sofia_income_noise.localMI) +
  tm_fill(col = "Pr.z....E.Ibvi..",
          breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
          palette = "-Blues",
          title = "local Moran's I p-values")+
  tm_borders(alpha = 0.5)

tmap_arrange(localMI.map, pvalue.map, asp = 1, ncol = 2)

### LISA cluster map
# Identifying the LISA clusters
quadrant <- vector(mode = "numeric", length = nrow(localMI))
signif <- 0.05 # significance level

# Assign quadrant values (HH, LL, HL, LH)
quadrant[sofia_income_noise.localMI$std_ses > 0 & sofia_income_noise.localMI$wsr_rate > 0] <- "High SES-High Mortality"      
quadrant[sofia_income_noise.localMI$std_ses < 0 & sofia_income_noise.localMI$wsr_rate < 0] <- "Low SES-Low Mortality"      
quadrant[sofia_income_noise.localMI$std_ses < 0 & sofia_income_noise.localMI$wsr_rate > 0] <- "Low SES-High Mortality"
quadrant[sofia_income_noise.localMI$std_ses > 0 & sofia_income_noise.localMI$wsr_rate < 0] <- "High SES-Low Mortality"
quadrant[localMI[,5]>signif] <- "NS"
sofia_income_noise.localMI$quadrant <- quadrant 

# Merge back to full polygon
sofia_income_noise.localMI$geometry <- NULL
sofia_income_noise <- merge(sofia_noise, sofia_income_noise.localMI[c("Poly_ID", "quadrant")], by = "Poly_ID", all.x = T)
sofia_income_noise <- sofia_income_noise %>% mutate(quadrant = replace_na(quadrant, "NS"))
sofia_income_noise$quadrant <- factor(sofia_income_noise$quadrant, levels = c("Low SES-High Mortality", "High SES-High Mortality", "High SES-Low Mortality", "Low SES-Low Mortality", "NS"))

# Make the plot
theme_set(theme_bw())
w8 <- ggplot(sofia_income_noise) + geom_sf(aes(fill = quadrant), size = 0.1) + 
  scale_fill_manual(name = "", values = c("High SES-High Mortality" = "#fdae61",  
                                          "Low SES-Low Mortality" = "#abd9e9", 
                                          "Low SES-High Mortality" = "#d7191c", 
                                          "High SES-Low Mortality" = "#2c7bb6", 
                                          "NS" = "#ffffff"), guide = guide_legend(nrow = 2, keywidth = 1.5, keyheight = 1.5)) +
  ggtitle("Noise mortality") + 
  theme(legend.position = "bottom", legend.text = element_text(size = 11), 
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        panel.grid = element_blank()) + 
  annotation_scale(location = "bl", width_hint = 0.2) +  # Add scale bar at bottom left
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering())

ggsave("Baseline_BOD/Sofia_spatial_analysis_Noise_mortality.png")


#-----------------------------------
### Green space mortality
sofia_green <- st_read("Baseline_BOD/Results/res_ndvi_nat_mort.geojson")
sofia_green <- st_transform(sofia_green, crs = 3035)
sofia_income_green <- merge(sofia_green, sofia_ses_poly[c("Poly_ID", "ses_index")], by = "Poly_ID")
sofia_income_green$ses_index <- -sofia_income_green$ses_index

# Plot data
qtm(sofia_income_green, "att_rate")
qtm(sofia_income_green, "ses_index")

# Spatial weights matrix
centroids <- st_centroid(st_geometry(sofia_income_green))
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
isolated_areas <- which(card(wm_d) == 0) # Identify areas with zero neighbors
sofia_income_green <- sofia_income_green[-isolated_areas, ] # Remove isolated areas from the data
centroids <- centroids[-isolated_areas, ]
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
# Row standardized weights matrix 
rswm_d <- nb2listw(wm_d, style = "W", zero.policy = TRUE)
summary(rswm_d)

### Analysis for green space
# Bivariate Moran I scatter plot
sofia_income_green$std_ses <- standardize(sofia_income_green$ses_index)
sofia_income_green$std_rate <- standardize(sofia_income_green$att_rate)
# Lag variable: 
sofia_income_green$wsr_rate <- lag.listw(rswm_d, sofia_income_green$std_rate)

# Plot: 
theme_set(theme_bw())
p5 <- ggplot(data = sofia_income_green, aes(x=std_ses, y=wsr_rate)) +
  geom_point(color = "grey", alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty =2) +
  ggtitle("SES and lagged NDVI mortality") + ylab("Standard lagged NDVI mortality") + xlab("Standard SES index")

ggsave("Baseline_BOD/Sofia_bivariateMoran_biplot_NDVI_mortality.png")

### Global Moran I (bivariate)
bivariate_moran <- moran_bv(sofia_income_green$std_ses, sofia_income_green$std_rate, rswm_d, nsim = 500)
print(bivariate_moran)
bivariate_moran$t0
summary(lm(sofia_income_green$wsr_rate ~ sofia_income_green$std_ses))

### Local Moran I (bivariate)
localMI <- localmoran_bv(sofia_income_green$std_ses, sofia_income_green$std_rate, rswm_d, nsim = 100)
head(localMI)

# Plot results
sofia_income_green.localMI <- cbind(sofia_income_green, localMI)

# Moran I
localMI.map <- tm_shape(sofia_income_green.localMI) +
  tm_fill(col = "Ibvi",
          style = "pretty",
          title = "local moran statistics")+
  tm_borders(alpha = 0.5)
# p values
pvalue.map <- tm_shape(sofia_income_green.localMI) +
  tm_fill(col = "Pr.z....E.Ibvi..",
          breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
          palette = "-Blues",
          title = "local Moran's I p-values")+
  tm_borders(alpha = 0.5)

tmap_arrange(localMI.map, pvalue.map, asp = 1, ncol = 2)

### LISA cluster map
# Identifying the LISA clusters
quadrant <- vector(mode = "numeric", length = nrow(localMI))
signif <- 0.05 # significance level

# Assign quadrant values (HH, LL, HL, LH)
quadrant[sofia_income_green.localMI$std_ses > 0 & sofia_income_green.localMI$wsr_rate > 0] <- "High SES-High Mortality"      
quadrant[sofia_income_green.localMI$std_ses < 0 & sofia_income_green.localMI$wsr_rate < 0] <- "Low SES-Low Mortality"      
quadrant[sofia_income_green.localMI$std_ses < 0 & sofia_income_green.localMI$wsr_rate > 0] <- "Low SES-High Mortality"
quadrant[sofia_income_green.localMI$std_ses > 0 & sofia_income_green.localMI$wsr_rate < 0] <- "High SES-Low Mortality"
quadrant[localMI[,5]>signif] <- "NS"
sofia_income_green.localMI$quadrant <- quadrant 

# Merge back to full polygon
sofia_income_green.localMI$geometry <- NULL
sofia_income_green <- merge(sofia_green, sofia_income_green.localMI[c("Poly_ID", "quadrant")], by = "Poly_ID", all.x = T)
sofia_income_green <- sofia_income_green %>% mutate(quadrant = replace_na(quadrant, "NS"))
sofia_income_green$quadrant <- factor(sofia_income_green$quadrant, levels = c("Low SES-High Mortality", "High SES-High Mortality", "High SES-Low Mortality", "Low SES-Low Mortality", "NS"))

# Make the plot
theme_set(theme_bw())
w10 <- ggplot(sofia_income_green) + geom_sf(aes(fill = quadrant), size = 0.1) + 
  scale_fill_manual(name = "", values = c("High SES-High Mortality" = "#fdae61",  
                                          "Low SES-Low Mortality" = "#abd9e9", 
                                          "Low SES-High Mortality" = "#d7191c", 
                                          "High SES-Low Mortality" = "#2c7bb6", 
                                          "NS" = "#ffffff"), guide = guide_legend(nrow = 2, keywidth = 1.5, keyheight = 1.5)) +
  ggtitle("Lack of green mortality") + 
  theme(legend.position = "bottom", legend.text = element_text(size = 11), 
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        panel.grid = element_blank()) + 
  annotation_scale(location = "bl", width_hint = 0.2) +  # Add scale bar at bottom left
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering())

ggsave("Baseline_BOD/Sofia_spatial_analysis_NDVI_mortality.png")


#-----------------------------------
### UHI mortality
sofia_uhi <- st_read("Baseline_BOD/Results/res_uhi_allcause_mort.geojson")
sofia_uhi <- st_transform(sofia_uhi, crs = 3035)
sofia_income_uhi <- merge(sofia_uhi, sofia_ses_poly[c("Poly_ID", "ses_index")], by = "Poly_ID")
sofia_income_uhi$ses_index <- -sofia_income_uhi$ses_index

# Plot data
qtm(sofia_income_uhi, "att_rate")
qtm(sofia_income_uhi, "ses_index")

# Spatial weights matrix
centroids <- st_centroid(st_geometry(sofia_income_uhi))
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
isolated_areas <- which(card(wm_d) == 0) # Identify areas with zero neighbors
sofia_income_uhi <- sofia_income_uhi[-isolated_areas, ] # Remove isolated areas from the data
centroids <- centroids[-isolated_areas, ]
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
# Row standardized weights matrix 
rswm_d <- nb2listw(wm_d, style = "W", zero.policy = TRUE)
summary(rswm_d)

### Analysis for UHI
# Bivariate Moran I scatter plot
sofia_income_uhi$std_ses <- standardize(sofia_income_uhi$ses_index)
sofia_income_uhi$std_rate <- standardize(sofia_income_uhi$att_rate)
# Lag variable: 
sofia_income_uhi$wsr_rate <- lag.listw(rswm_d,sofia_income_uhi$std_rate)

# Plot: 
theme_set(theme_bw())
p6 <- ggplot(data = sofia_income_uhi, aes(x=std_ses,y=wsr_rate)) +
  geom_point(color = "grey", alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty =2) +
  ggtitle("SES and lagged UHI mortality") + ylab("Standard lagged UHI mortality") + xlab("Standard SES index")

### Global Moran I (bivariate)
bivariate_moran <- moran_bv(sofia_income_uhi$std_ses, sofia_income_uhi$std_rate, rswm_d, nsim = 500)
print(bivariate_moran)
bivariate_moran$t0
summary(lm(sofia_income_uhi$wsr_rate ~ sofia_income_uhi$std_ses))

### Local Moran I (bivariate)
localMI <- localmoran_bv(sofia_income_uhi$std_ses, sofia_income_uhi$std_rate, rswm_d, nsim = 100)
head(localMI)

# Plot results
sofia_income_uhi.localMI <- cbind(sofia_income_uhi, localMI)

# Moran I
localMI.map <- tm_shape(sofia_income_uhi.localMI) +
  tm_fill(col = "Ibvi",
          style = "pretty",
          title = "local moran statistics")+
  tm_borders(alpha = 0.5)
# p values
pvalue.map <- tm_shape(sofia_income_uhi.localMI) +
  tm_fill(col = "Pr.z....E.Ibvi..",
          breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
          palette = "-Blues",
          title = "local Moran's I p-values")+
  tm_borders(alpha = 0.5)

tmap_arrange(localMI.map, pvalue.map, asp = 1, ncol = 2)

### LISA cluster map
# Identifying the LISA clusters
quadrant <- vector(mode = "numeric", length = nrow(localMI))
signif <- 0.05 # significance level

# Assign quadrant values (HH, LL, HL, LH)
quadrant[sofia_income_uhi.localMI$std_ses > 0 & sofia_income_uhi.localMI$wsr_rate > 0] <- "High SES-High Mortality"      
quadrant[sofia_income_uhi.localMI$std_ses < 0 & sofia_income_uhi.localMI$wsr_rate < 0] <- "Low SES-Low Mortality"      
quadrant[sofia_income_uhi.localMI$std_ses < 0 & sofia_income_uhi.localMI$wsr_rate > 0] <- "Low SES-High Mortality"
quadrant[sofia_income_uhi.localMI$std_ses > 0 & sofia_income_uhi.localMI$wsr_rate < 0] <- "High SES-Low Mortality"
quadrant[localMI[,5]>signif] <- "NS"
sofia_income_uhi.localMI$quadrant <- quadrant 

# Merge back to full polygon
sofia_income_uhi.localMI$geometry <- NULL
sofia_income_uhi <- merge(sofia_uhi, sofia_income_uhi.localMI[c("Poly_ID", "quadrant")], by = "Poly_ID", all.x = T)
sofia_income_uhi <- sofia_income_uhi %>% mutate(quadrant = replace_na(quadrant, "NS"))
sofia_income_uhi$quadrant <- factor(sofia_income_uhi$quadrant, levels = c("Low SES-High Mortality", "High SES-High Mortality", "High SES-Low Mortality", "Low SES-Low Mortality", "NS"))

# Make the plot
theme_set(theme_bw())
w12 <- ggplot(sofia_income_uhi) + geom_sf(aes(fill = quadrant), size = 0.1) + 
  scale_fill_manual(name = "", values = c("High SES-High Mortality" = "#fdae61",  
                                          "Low SES-Low Mortality" = "#abd9e9", 
                                          "Low SES-High Mortality" = "#d7191c", 
                                          "High SES-Low Mortality" = "#2c7bb6", 
                                          "NS" = "#ffffff"), guide = guide_legend(nrow = 2, keywidth = 1.5, keyheight = 1.5)) +
  ggtitle("UHI mortality") + 
  theme(legend.position = "bottom", legend.text = element_text(size = 11), 
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        panel.grid = element_blank()) + 
  annotation_scale(location = "bl", width_hint = 0.2) +  # Add scale bar at bottom left
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering())

ggsave("Baseline_BOD/Sofia_spatial_analysis_UHI_mortality.png")


### Plot with all mortality results
library(ggpubr)
ggarrange(w2, w4, w6, w8, w10, w12, ncol = 3, nrow = 2)
ggsave("Baseline_BOD/Sofia_spatial_analysis_polygons_mortality.png")

# Bivariate Moran I biplot
ggarrange(p1, p2, p3, p4, p5, p6,  ncol = 2, nrow = 3)
ggsave("Baseline_BOD/Sofia_bivariateMoran_biplot_mortality.png")


#-----------------------------------
### PM2.5 IHD incidence
sofia_airp <- st_read("Baseline_BOD/Results/res_pm25_SES_ihd.geojson")
sofia_airp <- st_transform(sofia_airp, crs = 3035)
sofia_income_airp <- merge(sofia_airp, sofia_ses_poly[c("Poly_ID", "ses_index")], by = "Poly_ID")
sofia_income_airp$ses_index <- -sofia_income_airp$ses_index

# Plot data
qtm(sofia_income_airp, "att_rate")
qtm(sofia_income_airp, "ses_index")

# Spatial weights matrix
centroids <- st_centroid(st_geometry(sofia_income_airp))
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
isolated_areas <- which(card(wm_d) == 0) # Identify areas with zero neighbors
sofia_income_airp <- sofia_income_airp[-isolated_areas, ] # Remove isolated areas from the data
centroids <- centroids[-isolated_areas, ]
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
# Row standardized weights matrix 
rswm_d <- nb2listw(wm_d, style = "W", zero.policy = TRUE)
summary(rswm_d)

### Analysis for PM2.5
# Bivariate Moran I scatter plot
sofia_income_airp$std_ses <- standardize(sofia_income_airp$ses_index)
sofia_income_airp$std_rate <- standardize(sofia_income_airp$att_rate)
# Lag variable: 
sofia_income_airp$wsr_rate <- lag.listw(rswm_d,sofia_income_airp$std_rate)

# Plot:
theme_set(theme_bw())
p1 <- ggplot(data = sofia_income_airp, aes(x=std_ses,y=wsr_rate)) +
  geom_point(color = "grey", alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty =2) +
  ggtitle(expression("SES and lagged PM"[2.5]*" IHD incidence")) + ylab(expression("Standard lagged PM"[2.5]*" IHD incidence")) + xlab("Standard SES index")

### Global Moran I (bivariate)
bivariate_moran <- moran_bv(sofia_income_airp$std_ses, sofia_income_airp$std_rate, rswm_d, nsim = 500)
print(bivariate_moran)
bivariate_moran$t0
summary(lm(sofia_income_airp$wsr_rate ~ sofia_income_airp$std_ses))

### Local Moran I (bivariate)
localMI <- localmoran_bv(sofia_income_airp$std_ses, sofia_income_airp$std_rate, rswm_d, nsim = 100)
head(localMI)

# Plot results
sofia_income_airp.localMI <- cbind(sofia_income_airp, localMI)

# Moran I
localMI.map <- tm_shape(sofia_income_airp.localMI) +
  tm_fill(col = "Ibvi",
          style = "pretty",
          title = "local moran statistics")+
  tm_borders(alpha = 0.5)
# p values
pvalue.map <- tm_shape(sofia_income_airp.localMI) +
  tm_fill(col = "Pr.z....E.Ibvi..",
          breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
          palette = "-Blues",
          title = "local Moran's I p-values")+
  tm_borders(alpha = 0.5)

tmap_arrange(localMI.map, pvalue.map, asp = 1, ncol = 2)

### LISA cluster map
# Identifying the LISA clusters
quadrant <- vector(mode = "numeric", length = nrow(localMI))
signif <- 0.05 # significance level

# Assign quadrant values (HH, LL, HL, LH)
quadrant[sofia_income_airp.localMI$std_ses > 0 & sofia_income_airp.localMI$wsr_rate > 0] <- "High SES-High Incidence"      
quadrant[sofia_income_airp.localMI$std_ses < 0 & sofia_income_airp.localMI$wsr_rate < 0] <- "Low SES-Low Incidence"      
quadrant[sofia_income_airp.localMI$std_ses < 0 & sofia_income_airp.localMI$wsr_rate > 0] <- "Low SES-High Incidence"
quadrant[sofia_income_airp.localMI$std_ses > 0 & sofia_income_airp.localMI$wsr_rate < 0] <- "High SES-Low Incidence"
quadrant[localMI[,5]>signif] <- "NS"
sofia_income_airp.localMI$quadrant <- quadrant 

# Merge back to full polygon
sofia_income_airp.localMI$geometry <- NULL
sofia_income_airp <- merge(sofia_airp, sofia_income_airp.localMI[c("Poly_ID", "quadrant")], by = "Poly_ID", all.x = T)
sofia_income_airp <- sofia_income_airp %>% mutate(quadrant = replace_na(quadrant, "NS"))
sofia_income_airp$quadrant <- factor(sofia_income_airp$quadrant, levels = c("Low SES-High Incidence", "High SES-High Incidence", "High SES-Low Incidence", "Low SES-Low Incidence", "NS"))

# Make the plot
theme_set(theme_bw())
w1 <- ggplot(sofia_income_airp) + geom_sf(aes(fill = quadrant), size = 0.1) + 
  scale_fill_manual(name = "", values = c("High SES-High Incidence" = "#fdae61",  
                                          "Low SES-Low Incidence" = "#abd9e9", 
                                          "Low SES-High Incidence" = "#d7191c", 
                                          "High SES-Low Incidence" = "#2c7bb6", 
                                          "NS" = "#ffffff"), guide = guide_legend(nrow = 2, keywidth = 1.5, keyheight = 1.5)) +
  ggtitle(expression("PM"[2.5]*" IHD incidence")) + 
  theme(legend.position = "bottom", legend.text = element_text(size = 11), 
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        panel.grid = element_blank()) + 
  annotation_scale(location = "bl", width_hint = 0.2) +  # Add scale bar at bottom left
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering())

ggsave("Baseline_BOD/Sofia_spatial_analysis_PM2.5_IHD.png")


#-----------------------------------
### PM2.5 COPD incidence
sofia_airp <- st_read("Baseline_BOD/Results/res_pm25_SES_copd.geojson")
sofia_airp <- st_transform(sofia_airp, crs = 3035)
sofia_income_airp <- merge(sofia_airp, sofia_ses_poly[c("Poly_ID", "ses_index")], by = "Poly_ID")
sofia_income_airp$ses_index <- -sofia_income_airp$ses_index

# Plot data
qtm(sofia_income_airp, "att_rate")
qtm(sofia_income_airp, "ses_index")

# Spatial weights matrix
centroids <- st_centroid(st_geometry(sofia_income_airp))
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
isolated_areas <- which(card(wm_d) == 0) # Identify areas with zero neighbors
sofia_income_airp <- sofia_income_airp[-isolated_areas, ] # Remove isolated areas from the data
centroids <- centroids[-isolated_areas, ]
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
# Row standardized weights matrix 
rswm_d <- nb2listw(wm_d, style = "W", zero.policy = TRUE)
summary(rswm_d)

### Analysis for NO2
# Bivariate Moran I scatter plot
sofia_income_airp$std_ses <- standardize(sofia_income_airp$ses_index)
sofia_income_airp$std_rate <- standardize(sofia_income_airp$att_rate)
# Lag variable: 
sofia_income_airp$wsr_rate <- lag.listw(rswm_d,sofia_income_airp$std_rate)

# Plot:
theme_set(theme_bw())
p2 <- ggplot(data = sofia_income_airp, aes(x=std_ses,y=wsr_rate)) +
  geom_point(color = "grey", alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty =2) +
  ggtitle(expression("SES and lagged PM"[2.5]*" COPD incidence")) + ylab(expression("Standard lagged PM"[2.5]*" COPD incidence")) + xlab("Standard SES index")

### Global Moran I (bivariate)
bivariate_moran <- moran_bv(sofia_income_airp$std_ses, sofia_income_airp$std_rate, rswm_d, nsim = 500)
print(bivariate_moran)
bivariate_moran$t0
summary(lm(sofia_income_airp$wsr_rate ~ sofia_income_airp$std_ses))

### Local Moran I (bivariate)
localMI <- localmoran_bv(sofia_income_airp$std_ses, sofia_income_airp$std_rate, rswm_d, nsim = 100)
head(localMI)

# Plot results
sofia_income_airp.localMI <- cbind(sofia_income_airp, localMI)

# Moran I
localMI.map <- tm_shape(sofia_income_airp.localMI) +
  tm_fill(col = "Ibvi",
          style = "pretty",
          title = "local moran statistics")+
  tm_borders(alpha = 0.5)
# p values
pvalue.map <- tm_shape(sofia_income_airp.localMI) +
  tm_fill(col = "Pr.z....E.Ibvi..",
          breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
          palette = "-Blues",
          title = "local Moran's I p-values")+
  tm_borders(alpha = 0.5)

tmap_arrange(localMI.map, pvalue.map, asp = 1, ncol = 2)

### LISA cluster map
# Identifying the LISA clusters
quadrant <- vector(mode = "numeric", length = nrow(localMI))
signif <- 0.05 # significance level

# Assign quadrant values (HH, LL, HL, LH)
quadrant[sofia_income_airp.localMI$std_ses > 0 & sofia_income_airp.localMI$wsr_rate > 0] <- "High SES-High Incidence"      
quadrant[sofia_income_airp.localMI$std_ses < 0 & sofia_income_airp.localMI$wsr_rate < 0] <- "Low SES-Low Incidence"      
quadrant[sofia_income_airp.localMI$std_ses < 0 & sofia_income_airp.localMI$wsr_rate > 0] <- "Low SES-High Incidence"
quadrant[sofia_income_airp.localMI$std_ses > 0 & sofia_income_airp.localMI$wsr_rate < 0] <- "High SES-Low Incidence"
quadrant[localMI[,5]>signif] <- "NS"
sofia_income_airp.localMI$quadrant <- quadrant 

# Merge back to full polygon
sofia_income_airp.localMI$geometry <- NULL
sofia_income_airp <- merge(sofia_airp, sofia_income_airp.localMI[c("Poly_ID", "quadrant")], by = "Poly_ID", all.x = T)
sofia_income_airp <- sofia_income_airp %>% mutate(quadrant = replace_na(quadrant, "NS"))
sofia_income_airp$quadrant <- factor(sofia_income_airp$quadrant, levels = c("Low SES-High Incidence", "High SES-High Incidence", "High SES-Low Incidence", "Low SES-Low Incidence", "NS"))

# Make the plot
theme_set(theme_bw())
w2 <- ggplot(sofia_income_airp) + geom_sf(aes(fill = quadrant), size = 0.1) + 
  scale_fill_manual(name = "", values = c("High SES-High Incidence" = "#fdae61",  
                                          "Low SES-Low Incidence" = "#abd9e9", 
                                          "Low SES-High Incidence" = "#d7191c", 
                                          "High SES-Low Incidence" = "#2c7bb6", 
                                          "NS" = "#ffffff"), guide = guide_legend(nrow = 2, keywidth = 1.5, keyheight = 1.5)) +
  ggtitle(expression("PM"[2.5]*" COPD incidence")) + 
  theme(legend.position = "bottom", legend.text = element_text(size = 11), 
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        panel.grid = element_blank()) + 
  annotation_scale(location = "bl", width_hint = 0.2) +  # Add scale bar at bottom left
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering())

ggsave("Baseline_BOD/Sofia_spatial_analysis_PM25_copd.png")


#-----------------------------------
### NO2 asthma in adults
sofia_airp <- st_read("Baseline_BOD/Results/res_no2_SES_asthma_adults.geojson")
sofia_airp <- st_transform(sofia_airp, crs = 3035)
sofia_income_airp <- merge(sofia_airp, sofia_ses_poly[c("Poly_ID", "ses_index")], by = "Poly_ID")
sofia_income_airp$ses_index <- -sofia_income_airp$ses_index

# Plot data
qtm(sofia_income_airp, "att_rate")
qtm(sofia_income_airp, "ses_index")

# Spatial weights matrix
centroids <- st_centroid(st_geometry(sofia_income_airp))
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
isolated_areas <- which(card(wm_d) == 0) # Identify areas with zero neighbors
sofia_income_airp <- sofia_income_airp[-isolated_areas, ] # Remove isolated areas from the data
centroids <- centroids[-isolated_areas, ]
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
# Row standardized weights matrix 
rswm_d <- nb2listw(wm_d, style = "W", zero.policy = TRUE)
summary(rswm_d)

### Analysis for NO2
# Bivariate Moran I scatter plot
sofia_income_airp$std_ses <- standardize(sofia_income_airp$ses_index)
sofia_income_airp$std_rate <- standardize(sofia_income_airp$att_rate)
# Lag variable: 
sofia_income_airp$wsr_rate <- lag.listw(rswm_d,sofia_income_airp$std_rate)

# Plot: 
theme_set(theme_bw())
p3 <- ggplot(data = sofia_income_airp, aes(x=std_ses,y=wsr_rate)) +
  geom_point(color = "grey", alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty =2) +
  ggtitle(expression("SES and lagged NO"[2]*" asthma incidence")) + ylab(expression("Standard lagged NO"[2]*" asthma incidence")) + xlab("Standard SES index")

### Global Moran I (bivariate)
bivariate_moran <- moran_bv(sofia_income_airp$std_ses, sofia_income_airp$std_rate, rswm_d, nsim = 500)
print(bivariate_moran)
bivariate_moran$t0
summary(lm(sofia_income_airp$wsr_rate ~ sofia_income_airp$std_ses))

### Local Moran I (bivariate)
localMI <- localmoran_bv(sofia_income_airp$std_ses, sofia_income_airp$std_rate, rswm_d, nsim = 100)
head(localMI)

# Plot results
sofia_income_airp.localMI <- cbind(sofia_income_airp, localMI)

# Moran I
localMI.map <- tm_shape(sofia_income_airp.localMI) +
  tm_fill(col = "Ibvi",
          style = "pretty",
          title = "local moran statistics")+
  tm_borders(alpha = 0.5)
# p values
pvalue.map <- tm_shape(sofia_income_airp.localMI) +
  tm_fill(col = "Pr.z....E.Ibvi..",
          breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
          palette = "-Blues",
          title = "local Moran's I p-values")+
  tm_borders(alpha = 0.5)

tmap_arrange(localMI.map, pvalue.map, asp = 1, ncol = 2)

### LISA cluster map
# Identifying the LISA clusters
quadrant <- vector(mode = "numeric", length = nrow(localMI))
signif <- 0.05 # significance level

# Assign quadrant values (HH, LL, HL, LH)
quadrant[sofia_income_airp.localMI$std_ses > 0 & sofia_income_airp.localMI$wsr_rate > 0] <- "High SES-High Incidence"      
quadrant[sofia_income_airp.localMI$std_ses < 0 & sofia_income_airp.localMI$wsr_rate < 0] <- "Low SES-Low Incidence"      
quadrant[sofia_income_airp.localMI$std_ses < 0 & sofia_income_airp.localMI$wsr_rate > 0] <- "Low SES-High Incidence"
quadrant[sofia_income_airp.localMI$std_ses > 0 & sofia_income_airp.localMI$wsr_rate < 0] <- "High SES-Low Incidence"
quadrant[localMI[,5]>signif] <- "NS"
sofia_income_airp.localMI$quadrant <- quadrant 

# Merge back to full polygon
sofia_income_airp.localMI$geometry <- NULL
sofia_income_airp <- merge(sofia_airp, sofia_income_airp.localMI[c("Poly_ID", "quadrant")], by = "Poly_ID", all.x = T)
sofia_income_airp <- sofia_income_airp %>% mutate(quadrant = replace_na(quadrant, "NS"))
sofia_income_airp$quadrant <- factor(sofia_income_airp$quadrant, levels = c("Low SES-High Incidence", "High SES-High Incidence", "High SES-Low Incidence", "Low SES-Low Incidence", "NS"))

# Make the plot
theme_set(theme_bw())
w3 <- ggplot(sofia_income_airp) + geom_sf(aes(fill = quadrant), size = 0.1) + 
  scale_fill_manual(name = "", values = c("High SES-High Incidence" = "#fdae61",  
                                          "Low SES-Low Incidence" = "#abd9e9", 
                                          "Low SES-High Incidence" = "#d7191c", 
                                          "High SES-Low Incidence" = "#2c7bb6", 
                                          "NS" = "#ffffff"), guide = guide_legend(nrow = 2, keywidth = 1.5, keyheight = 1.5)) +
  ggtitle(expression("NO"[2]*" asthma incidence in adults")) + 
  theme(legend.position = "bottom", legend.text = element_text(size = 11), 
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        panel.grid = element_blank()) + 
  annotation_scale(location = "bl", width_hint = 0.2) +  # Add scale bar at bottom left
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering())

ggsave("Baseline_BOD/Sofia_spatial_analysis_NO2_asthma_adults.png")


#-----------------------------------
### Noise IHD incidence
sofia_noise <- st_read("Baseline_BOD/Results/res_noise_SES_ihd.geojson")
sofia_noise <- st_transform(sofia_noise, crs = 3035)
sofia_income_noise <- merge(sofia_noise, sofia_ses_poly[c("Poly_ID", "ses_index")], by = "Poly_ID")
sofia_income_noise$ses_index <- -sofia_income_noise$ses_index

# Plot data
qtm(sofia_income_noise, "att_rate")
qtm(sofia_income_noise, "ses_index")

# Spatial weights matrix
centroids <- st_centroid(st_geometry(sofia_income_noise))
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
isolated_areas <- which(card(wm_d) == 0) # Identify areas with zero neighbors
sofia_income_noise <- sofia_income_noise[-isolated_areas, ] # Remove isolated areas from the data
centroids <- centroids[-isolated_areas, ]
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
# Row standardized weights matrix 
rswm_d <- nb2listw(wm_d, style = "W", zero.policy = TRUE)
summary(rswm_d)

### Analysis for noise
# Bivariate Moran I scatter plot
sofia_income_noise$std_ses <- standardize(sofia_income_noise$ses_index)
sofia_income_noise$std_rate <- standardize(sofia_income_noise$att_rate)
# Lag variable: 
sofia_income_noise$wsr_rate <- lag.listw(rswm_d,sofia_income_noise$std_rate)

# Plot: 
theme_set(theme_bw())
p4 <- ggplot(data = sofia_income_noise, aes(x=std_ses,y=wsr_rate)) +
  geom_point(color = "grey", alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty =2) +
  ggtitle("SES and lagged noise IHD incidence") + ylab("Standard lagged noise IHD incidence") + xlab("Standard SES index")

### Global Moran I (bivariate)
bivariate_moran <- moran_bv(sofia_income_noise$std_ses, sofia_income_noise$std_rate, rswm_d, nsim = 500)
print(bivariate_moran)
bivariate_moran$t0
summary(lm(sofia_income_noise$wsr_rate ~ sofia_income_noise$std_ses))

### Local Moran I (bivariate)
localMI <- localmoran_bv(sofia_income_noise$std_ses, sofia_income_noise$std_rate, rswm_d, nsim = 100)
head(localMI)

# Plot results
sofia_income_noise.localMI <- cbind(sofia_income_noise, localMI)

# Moran I
localMI.map <- tm_shape(sofia_income_noise.localMI) +
  tm_fill(col = "Ibvi",
          style = "pretty",
          title = "local moran statistics")+
  tm_borders(alpha = 0.5)
# p values
pvalue.map <- tm_shape(sofia_income_noise.localMI) +
  tm_fill(col = "Pr.z....E.Ibvi..",
          breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
          palette = "-Blues",
          title = "local Moran's I p-values")+
  tm_borders(alpha = 0.5)

tmap_arrange(localMI.map, pvalue.map, asp = 1, ncol = 2)

### LISA cluster map
# Identifying the LISA clusters
quadrant <- vector(mode = "numeric", length = nrow(localMI))
signif <- 0.05 # significance level

# Assign quadrant values (HH, LL, HL, LH)
quadrant[sofia_income_noise.localMI$std_ses > 0 & sofia_income_noise.localMI$wsr_rate > 0] <- "High SES-High Incidence"      
quadrant[sofia_income_noise.localMI$std_ses < 0 & sofia_income_noise.localMI$wsr_rate < 0] <- "Low SES-Low Incidence"      
quadrant[sofia_income_noise.localMI$std_ses < 0 & sofia_income_noise.localMI$wsr_rate > 0] <- "Low SES-High Incidence"
quadrant[sofia_income_noise.localMI$std_ses > 0 & sofia_income_noise.localMI$wsr_rate < 0] <- "High SES-Low Incidence"
quadrant[localMI[,5]>signif] <- "NS"
sofia_income_noise.localMI$quadrant <- quadrant 

# Merge back to full polygon
sofia_income_noise.localMI$geometry <- NULL
sofia_income_noise <- merge(sofia_noise, sofia_income_noise.localMI[c("Poly_ID", "quadrant")], by = "Poly_ID", all.x = T)
sofia_income_noise <- sofia_income_noise %>% mutate(quadrant = replace_na(quadrant, "NS"))
sofia_income_noise$quadrant <- factor(sofia_income_noise$quadrant, levels = c("Low SES-High Incidence", "High SES-High Incidence", "High SES-Low Incidence", "Low SES-Low Incidence", "NS"))

# Make the plot
theme_set(theme_bw())
w4 <- ggplot(sofia_income_noise) + geom_sf(aes(fill = quadrant), size = 0.1) + 
  scale_fill_manual(name = "", values = c("High SES-High Incidence" = "#fdae61",  
                                          "Low SES-Low Incidence" = "#abd9e9", 
                                          "Low SES-High Incidence" = "#d7191c", 
                                          "High SES-Low Incidence" = "#2c7bb6", 
                                          "NS" = "#ffffff"), guide = guide_legend(nrow = 2, keywidth = 1.5, keyheight = 1.5)) +
  ggtitle("Noise IHD incidence") + 
  theme(legend.position = "bottom", legend.text = element_text(size = 11), 
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        panel.grid = element_blank()) + 
  annotation_scale(location = "bl", width_hint = 0.2) +  # Add scale bar at bottom left
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering())

ggsave("Baseline_BOD/Sofia_spatial_analysis_Noise_IHD.png")


#-----------------------------------
### Noise diabetes type II incidence
sofia_noise <- st_read("Baseline_BOD/Results/res_noise_SES_diabetes.geojson")
sofia_noise <- st_transform(sofia_noise, crs = 3035)
sofia_income_noise <- merge(sofia_noise, sofia_ses_poly[c("Poly_ID", "ses_index")], by = "Poly_ID")
sofia_income_noise$ses_index <- -sofia_income_noise$ses_index

# Plot data
qtm(sofia_income_noise, "att_rate")
qtm(sofia_income_noise, "ses_index")

# Spatial weights matrix
centroids <- st_centroid(st_geometry(sofia_income_noise))
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
isolated_areas <- which(card(wm_d) == 0) # Identify areas with zero neighbors
sofia_income_noise <- sofia_income_noise[-isolated_areas, ] # Remove isolated areas from the data
centroids <- centroids[-isolated_areas, ]
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
# Row standardized weights matrix 
rswm_d <- nb2listw(wm_d, style = "W", zero.policy = TRUE)
summary(rswm_d)

### Analysis for noise
# Bivariate Moran I scatter plot
sofia_income_noise$std_ses <- standardize(sofia_income_noise$ses_index)
sofia_income_noise$std_rate <- standardize(sofia_income_noise$att_rate)
# Lag variable: 
sofia_income_noise$wsr_rate <- lag.listw(rswm_d,sofia_income_noise$std_rate)

# Plot: 
theme_set(theme_bw())
p5 <- ggplot(data = sofia_income_noise, aes(x=std_ses,y=wsr_rate)) +
  geom_point(color = "grey", alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty =2) +
  ggtitle("SES and lagged noise diabetes (type II) incidence") + ylab("Standard lagged noise diabetes incidence") + xlab("Standard SES index")

### Global Moran I (bivariate)
bivariate_moran <- moran_bv(sofia_income_noise$std_ses, sofia_income_noise$std_rate, rswm_d, nsim = 500)
print(bivariate_moran)
bivariate_moran$t0
summary(lm(sofia_income_noise$wsr_rate ~ sofia_income_noise$std_ses))

### Local Moran I (bivariate)
localMI <- localmoran_bv(sofia_income_noise$std_ses, sofia_income_noise$std_rate, rswm_d, nsim = 100)
head(localMI)

# Plot results
sofia_income_noise.localMI <- cbind(sofia_income_noise, localMI)

# Moran I
localMI.map <- tm_shape(sofia_income_noise.localMI) +
  tm_fill(col = "Ibvi",
          style = "pretty",
          title = "local moran statistics")+
  tm_borders(alpha = 0.5)
# p values
pvalue.map <- tm_shape(sofia_income_noise.localMI) +
  tm_fill(col = "Pr.z....E.Ibvi..",
          breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
          palette = "-Blues",
          title = "local Moran's I p-values")+
  tm_borders(alpha = 0.5)

tmap_arrange(localMI.map, pvalue.map, asp = 1, ncol = 2)

### LISA cluster map
# Identifying the LISA clusters
quadrant <- vector(mode = "numeric", length = nrow(localMI))
signif <- 0.05 # significance level

# Assign quadrant values (HH, LL, HL, LH)
quadrant[sofia_income_noise.localMI$std_ses > 0 & sofia_income_noise.localMI$wsr_rate > 0] <- "High SES-High Incidence"      
quadrant[sofia_income_noise.localMI$std_ses < 0 & sofia_income_noise.localMI$wsr_rate < 0] <- "Low SES-Low Incidence"      
quadrant[sofia_income_noise.localMI$std_ses < 0 & sofia_income_noise.localMI$wsr_rate > 0] <- "Low SES-High Incidence"
quadrant[sofia_income_noise.localMI$std_ses > 0 & sofia_income_noise.localMI$wsr_rate < 0] <- "High SES-Low Incidence"
quadrant[localMI[,5]>signif] <- "NS"
sofia_income_noise.localMI$quadrant <- quadrant 

# Merge back to full polygon
sofia_income_noise.localMI$geometry <- NULL
sofia_income_noise <- merge(sofia_noise, sofia_income_noise.localMI[c("Poly_ID", "quadrant")], by = "Poly_ID", all.x = T)
sofia_income_noise <- sofia_income_noise %>% mutate(quadrant = replace_na(quadrant, "NS"))
sofia_income_noise$quadrant <- factor(sofia_income_noise$quadrant, levels = c("Low SES-High Incidence", "High SES-High Incidence", "High SES-Low Incidence", "Low SES-Low Incidence", "NS"))

# Make the plot
theme_set(theme_bw())
w5 <- ggplot(sofia_income_noise) + geom_sf(aes(fill = quadrant), size = 0.1) + 
  scale_fill_manual(name = "", values = c("High SES-High Incidence" = "#fdae61",  
                                          "Low SES-Low Incidence" = "#abd9e9", 
                                          "Low SES-High Incidence" = "#d7191c", 
                                          "High SES-Low Incidence" = "#2c7bb6", 
                                          "NS" = "#ffffff"), guide = guide_legend(nrow = 2, keywidth = 1.5, keyheight = 1.5)) +
  ggtitle("Noise diabetes (type II) incidence") + 
  theme(legend.position = "bottom", legend.text = element_text(size = 11), 
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        panel.grid = element_blank()) + 
  annotation_scale(location = "bl", width_hint = 0.2) +  # Add scale bar at bottom left
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering())

ggsave("Baseline_BOD/Sofia_spatial_analysis_Noise_diabetes.png")


#-----------------------------------
### Green space stroke incidence
sofia_green <- st_read("Baseline_BOD/Results/res_ndvi_SES_stroke.geojson")
sofia_green <- st_transform(sofia_green, crs = 3035)
sofia_income_green <- merge(sofia_green, sofia_ses_poly[c("Poly_ID", "ses_index")], by = "Poly_ID")
sofia_income_green$ses_index <- -sofia_income_green$ses_index

# Plot data
qtm(sofia_income_green, "att_rate")
qtm(sofia_income_green, "ses_index")

# Spatial weights matrix
centroids <- st_centroid(st_geometry(sofia_income_green))
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
isolated_areas <- which(card(wm_d) == 0) # Identify areas with zero neighbors
sofia_income_green <- sofia_income_green[-isolated_areas, ] # Remove isolated areas from the data
centroids <- centroids[-isolated_areas, ]
wm_d <- dnearneigh(st_coordinates(centroids), 0, 500)
summary(wm_d)
# Row standardized weights matrix 
rswm_d <- nb2listw(wm_d, style = "W", zero.policy = TRUE)
summary(rswm_d)

### Analysis for green space
# Bivariate Moran I scatter plot
sofia_income_green$std_ses <- standardize(sofia_income_green$ses_index)
sofia_income_green$std_rate <- standardize(sofia_income_green$att_rate)
# Lag variable: 
sofia_income_green$wsr_rate <- lag.listw(rswm_d, sofia_income_green$std_rate)

# Plot: 
theme_set(theme_bw())
p6 <- ggplot(data = sofia_income_green, aes(x=std_ses, y=wsr_rate)) +
  geom_point(color = "grey", alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty =2) +
  ggtitle("SES and lagged NDVI stroke incidence") + ylab("Standard lagged NDVI stroke incidence") + xlab("Standard SES index")

### Global Moran I (bivariate)
bivariate_moran <- moran_bv(sofia_income_green$std_ses, sofia_income_green$std_rate, rswm_d, nsim = 500)
print(bivariate_moran)
bivariate_moran$t0
summary(lm(sofia_income_green$wsr_rate ~ sofia_income_green$std_ses))

### Local Moran I (bivariate)
localMI <- localmoran_bv(sofia_income_green$std_ses, sofia_income_green$std_rate, rswm_d, nsim = 100)
head(localMI)

# Plot results
sofia_income_green.localMI <- cbind(sofia_income_green, localMI)

# Moran I
localMI.map <- tm_shape(sofia_income_green.localMI) +
  tm_fill(col = "Ibvi",
          style = "pretty",
          title = "local moran statistics")+
  tm_borders(alpha = 0.5)
# p values
pvalue.map <- tm_shape(sofia_income_green.localMI) +
  tm_fill(col = "Pr.z....E.Ibvi..",
          breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
          palette = "-Blues",
          title = "local Moran's I p-values")+
  tm_borders(alpha = 0.5)

tmap_arrange(localMI.map, pvalue.map, asp = 1, ncol = 2)

### LISA cluster map
# Identifying the LISA clusters
quadrant <- vector(mode = "numeric", length = nrow(localMI))
signif <- 0.05 # significance level

# Assign quadrant values (HH, LL, HL, LH)
quadrant[sofia_income_green.localMI$std_ses > 0 & sofia_income_green.localMI$wsr_rate > 0] <- "High SES-High Incidence"      
quadrant[sofia_income_green.localMI$std_ses < 0 & sofia_income_green.localMI$wsr_rate < 0] <- "Low SES-Low Incidence"      
quadrant[sofia_income_green.localMI$std_ses < 0 & sofia_income_green.localMI$wsr_rate > 0] <- "Low SES-High Incidence"
quadrant[sofia_income_green.localMI$std_ses > 0 & sofia_income_green.localMI$wsr_rate < 0] <- "High SES-Low Incidence"
quadrant[localMI[,5]>signif] <- "NS"
sofia_income_green.localMI$quadrant <- quadrant 

# Merge back to full polygon
sofia_income_green.localMI$geometry <- NULL
sofia_income_green <- merge(sofia_green, sofia_income_green.localMI[c("Poly_ID", "quadrant")], by = "Poly_ID", all.x = T)
sofia_income_green <- sofia_income_green %>% mutate(quadrant = replace_na(quadrant, "NS"))
sofia_income_green$quadrant <- factor(sofia_income_green$quadrant, levels = c("Low SES-High Incidence", "High SES-High Incidence", "High SES-Low Incidence", "Low SES-Low Incidence", "NS"))

# Make the plot
theme_set(theme_bw())
w6 <- ggplot(sofia_income_green) + geom_sf(aes(fill = quadrant), size = 0.1) + 
  scale_fill_manual(name = "", values = c("High SES-High Incidence" = "#fdae61",  
                                          "Low SES-Low Incidence" = "#abd9e9", 
                                          "Low SES-High Incidence" = "#d7191c", 
                                          "High SES-Low Incidence" = "#2c7bb6", 
                                          "NS" = "#ffffff"), guide = guide_legend(nrow = 2, keywidth = 1.5, keyheight = 1.5)) +
  ggtitle("Lack of green stroke incidence") + 
  theme(legend.position = "bottom", legend.text = element_text(size = 11), 
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        panel.grid = element_blank()) + 
  annotation_scale(location = "bl", width_hint = 0.2) +  # Add scale bar at bottom left
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering())

ggsave("Baseline_BOD/Sofia_spatial_analysis_NDVI_stroke.png")

### Plot with all morbidity results
library(ggpubr)
ggarrange(w1, w2, w3, w4, w5, w6, ncol = 3, nrow = 2)
ggsave("Baseline_BOD/Sofia_spatial_analysis_polygons_morbidity_SESadj.png")

# Bivariate Moran I biplot
ggarrange(p1, p2, p3, p4, p5, p6,  ncol = 2, nrow = 3)
ggsave("Baseline_BOD/Sofia_bivariateMoran_biplot_morbidity_SESadj.png")


# # Join with the corresponding urban unit
# sofia_ct <- st_read("Baseline_BOD/neighborhoods_20200616.geojson")
# sofia_ct <- st_transform(sofia_ct, st_crs(sofia_income_airp))
# sofia_income_airp$area <- st_area(sofia_income_airp)
# sofia_income_airp <- st_intersection(sofia_income_airp, sofia_ct)
# sofia_income_airp$area2 <- st_area(sofia_income_airp)
# sofia_income_airp$perarea <- as.numeric((sofia_income_airp$area2 / sofia_income_airp$area) * 100)
# sofia_income_airp[duplicated(sofia_income_airp$Poly_ID),]$Poly_ID
# sofia_income_airp <- sofia_income_airp %>% group_by(Poly_ID) %>% filter(perarea==max(perarea))
# sofia_income_airp$geometry <- NULL
# sofia_income_airp$area <- NULL
# sofia_income_airp$area2 <- NULL
# sofia_income_airp$perarea <- NULL
# 
# # Summary by urban unit
# sofia_income_airp <- sofia_income_airp %>% group_by(id) %>% summarise(pm25 = mean(pm25, na.rm = T), no2 = mean(no2, na.rm = T), ses_index = mean(ses_index, na.rm = T))
# sofia_income_airp <- merge(sofia_ct, sofia_income_airp, by = "id")

# wm_q <- poly2nb(sofia_income_airp, queen = TRUE)
#summary(wm_q)
# Row standardized weights matrix 
# rswm_q <- nb2listw(wm_q, zero.policy = TRUE)
# summary(rswm_q)

# w1 <- tmap_grob(w1)
# w2 <- tmap_grob(w2)
# w3 <- tmap_grob(w3)
# w4 <- tmap_grob(w4)
# w5 <- tmap_grob(w5)
# w6 <- tmap_grob(w6)
# w7 <- tmap_grob(w7)
# w8 <- tmap_grob(w8)
# w9 <- tmap_grob(w9)
# w10 <- tmap_grob(w10)
