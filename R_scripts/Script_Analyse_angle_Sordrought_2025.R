

library(tidyr)
library(dplyr)
library(ggplot2)
library(agricolae)
library(multcompView)
library(statgenSTA)
library(dplyr)
library(readxl)


# setwd("C:/Users/2025an002/Desktop/Sordrought_2025/Aida_GWAS/Analyse_root_angle/data")
# list.files() 
data <- read_excel("data/pheno/data_angle_sordrought_2025.xlsx", sheet = 2)

geno_id <- unique(data$LIGNEE)
cross_id <- substr(x = data$LIGNEE, 1, 4)
table(cross_id)

# option add cross
# data$cross <- cross_id

str(data)

data$Root_angle <- as.numeric(data$Root_angle) 
#data <- na.omit(data)


#L? on obtient pour chaque g?notype, les donn?es ab?rentes pour l'ensemble des r?p?titions
data_angle <- data %>% 
  ggplot(aes(x = LIGNEE, y = Root_angle)) +
  geom_boxplot(aes(x = LIGNEE, y = Root_angle), outlier.shape = 8, outlier.color = "red", outlier.size = 3) +  
  geom_jitter(width = 0.2, alpha = 0.7, color = "blue") + 
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
data_angle


####### Fonction remove outliers
remove_outliers <- function(data, column) {
  Q1 <- quantile(data[[column]], 0.25, na.rm = TRUE)
  Q3 <- quantile(data[[column]], 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  
  data_filtered <- data[data[[column]] >= lower_bound & data[[column]] <= upper_bound, ]
  return(data_filtered)
}

# Loop through each Genotype and remove outliers
Genotypes <- unique(data$LIGNEE)
data_angle_cleaned <- data.frame()  # Empty dataframe to store results

for (gen in Genotypes) {
  subset_data <- data %>% filter(LIGNEE == gen)  # Subset for each Genotype
  cleaned_data <- remove_outliers(subset_data, "Root_angle")  # Apply outlier removal
  data_angle_cleaned <- rbind(data_angle_cleaned, cleaned_data)  # Append cleaned data
}

print(data_angle_cleaned)

#write.table(data_Root_angle_clean, "data_Root_angle_clean.txt", sep = ";", row.names = FALSE, col.names = TRUE)
#Ici les outliers sont ?limin?s

data_angle_wo_Out <- data_angle_cleaned %>% 
  ggplot(aes(x = LIGNEE, y = Root_angle)) +
  geom_boxplot(aes(x = LIGNEE, y = Root_angle), outlier.shape = 8, outlier.color = "red", outlier.size = 3) +  
  geom_jitter(width = 0.2, alpha = 0.7, color = "blue") + 
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
data_angle_wo_Out


cat("Avant :", nrow(data), "\n")
cat("Après :", nrow(data_angle_cleaned), "\n")
cat("Supprimées :", nrow(data) - nrow(data_angle_cleaned), "\n")


#Avant : 2112 
#Après : 1988 
#Supprimées : 124 

##### Compute mean
# on calcule la moyenne des trois plants pour chaque g?notype. Chaque point repr?sente une r?p?tition

data_angle_wo_Out_means <- data_angle_cleaned %>%
  group_by(Y, X, X_Y, Block, Rep, Plot, LIGNEE) %>%
  summarise(Root_angle = mean(Root_angle), .groups = "drop")



data_angle_c <- data_angle_wo_Out_means %>% 
  ggplot(aes(x = LIGNEE, y = Root_angle)) +
  geom_boxplot(aes(x = LIGNEE, y = Root_angle), outlier.shape = 8, outlier.color = "red", outlier.size = 3) +  
  geom_jitter(width = 0.2, alpha = 0.7, color = "blue") + 
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
data_angle_c

########################################################
############# Correction StatGen#######################
################ On raw data###########################
########################################################

data_root_Mean_angle <- data_angle_wo_Out_means

data_root_Mean_angle$rowId <- as.factor(data_root_Mean_angle$Y)
data_root_Mean_angle$colId <- as.factor(data_root_Mean_angle$X)
data_root_Mean_angle$LIGNEE <- as.character(data_root_Mean_angle$LIGNEE)
data_root_Mean_angle$Root_angle <- as.numeric(data_root_Mean_angle$Root_angle)

data_root_Mean_angle <- data_root_Mean_angle %>%
  filter(
    is.finite(Root_angle),
    !is.na(X),
    !is.na(Y)
  )

data_root_Mean_angle_R <- data_root_Mean_angle %>%
  dplyr::select(LIGNEE, X, Y, rowId, colId, Block, Rep, Plot, Root_angle)



data_root_Mean_angle_R <- createTD(data = data_root_Mean_angle_R,
                                   genotype = "LIGNEE",
                                   rowCoord = "Y",
                                   colCoord = "X",
                                   rowId = "rowId",
                                   colId = "colId",
                                   subBlock = "Block",
                                   repId = "Rep",
                                   plotId = "Plot")

### Geno random
spaMod_data_angle_R <- fitTD(TD = data_root_Mean_angle_R,
                             design = "res.ibd",
                             traits = "Root_angle",
                             what = "random")

summary(spaMod_data_angle_R)

## Create spatial plots of the results.
plot(spaMod_data_angle_R, plotType = "spatial", spaTrend = ("percentage"))

## Detect outliers in the standardized residuals of the fitted model.
spats_outliers <- outlierSTA(STA = spaMod_data_angle_R, traits = "Root_angle", what = "random")
spats_outliers

## Extract all available statistics from the fitted model.
spats_extr <- extractSTA(spaMod_data_angle_R, what = "heritability")
spats_extr


#SUPPRESSION DES OUTLIERS
data_root_Mean_angle_R <- data_root_Mean_angle %>%
  dplyr::select(LIGNEE, X, Y, rowId, colId, Block, Rep, Plot, Root_angle)



data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="179")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="562")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="412")


data_root_Mean_angle_R <- createTD(data = data_root_Mean_angle_R,
                                   genotype = "LIGNEE",
                                   rowCoord = "Y",
                                   colCoord = "X",
                                   rowId = "rowId",
                                   colId = "colId",
                                   subBlock = "Block",
                                   repId = "Rep",
                                   plotId = "Plot")

### Geno random
spaMod_data_angle_R <- fitTD(TD = data_root_Mean_angle_R,
                             design = "res.ibd",
                             traits = "Root_angle",
                             what = "random")

summary(spaMod_data_angle_R)

## Create spatial plots of the results.
plot(spaMod_data_angle_R, plotType = "spatial", spaTrend = ("percentage"))

## Detect outliers in the standardized residuals of the fitted model.
spats_outliers <- outlierSTA(STA = spaMod_data_angle_R, traits = "Root_angle", what = "random")
spats_outliers

## Extract all available statistics from the fitted model.
spats_extr <- extractSTA(spaMod_data_angle_R, what = "heritability")
spats_extr


### Geno fixed
data_root_Mean_angle_F <- data_root_Mean_angle %>%
  select(LIGNEE, X, Y, rowId, colId, Block, Rep, Plot, Root_angle)


data_root_Mean_angle_F <- createTD(data = data_root_Mean_angle_F,
                                   genotype = "LIGNEE",
                                   rowCoord = "Y",
                                   colCoord = "X",
                                   rowId = "rowId",
                                   colId = "colId",
                                   subBlock = "Block",
                                   repId = "Rep",
                                   plotId = "Plot")

spaMod_Biom_F <- fitTD(TD = data_root_Mean_angle_F,
                       design = "res.ibd",
                       traits = "Root_angle",
                       what = "fixed")

summary(spaMod_Biom_F)

## Create spatial plots of the results.
plot(spaMod_Biom_F, plotType = "spatial", spaTrend = ("percentage"))

## Detect outliers in the standardized residuals of the fitted model.
spats_outliers <- outlierSTA(STA = spaMod_Biom_F, traits = "Root_angle", what = "fixed")
spats_outliers

#extr_lme4 <- extractSTA(spaMod, what=c("sMEANs"),  restoreColNames = TRUE, keep=c("repId","plotId","subBlock","rowCoord","colCoord"))
spats_TDGxEf <- STAtoTD(STA = spaMod_Biom_F, what = c("BLUEs", "seBLUEs"))
spats_TDGxEf




#Suppression des outliers

data_root_Mean_angle_F <- data_root_Mean_angle %>%
  select(LIGNEE, X, Y, rowId, colId, Block, Rep, Plot, Root_angle)

data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="276")
data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="412")


data_root_Mean_angle_F <- createTD(data = data_root_Mean_angle_F,
                                   genotype = "LIGNEE",
                                   rowCoord = "Y",
                                   colCoord = "X",
                                   rowId = "rowId",
                                   colId = "colId",
                                   subBlock = "Block",
                                   repId = "Rep",
                                   plotId = "Plot")

spaMod_Biom_F <- fitTD(TD = data_root_Mean_angle_F,
                       design = "res.ibd",
                       traits = "Root_angle",
                       what = "fixed")

summary(spaMod_Biom_F)

## Create spatial plots of the results.
plot(spaMod_Biom_F, plotType = "spatial", spaTrend = ("percentage"))

## Detect outliers in the standardized residuals of the fitted model.
spats_outliers <- outlierSTA(STA = spaMod_Biom_F, traits = "Root_angle", what = "fixed")
spats_outliers
#extr_lme4 <- extractSTA(spaMod, what=c("sMEANs"),  restoreColNames = TRUE, keep=c("repId","plotId","subBlock","rowCoord","colCoord"))
spats_TDGxEf <- STAtoTD(STA = spaMod_Biom_F, what = c("BLUEs", "seBLUEs"))
spats_TDGxEf

#Write blups 
write.table(spats_TDGxEf, "BLUEs_Mean_angle_geno.txt", sep = ";", row.names = FALSE, col.names = TRUE)
write.csv(spats_TDGxEf, "BLUEs_Mean_angle_geno.csv")

write.csv(spats_TDGxEf, "data/pheno/BLUEs_Mean_angle_geno.csv")


#Pour enrégistrer les BLUES
str(spats_TDGxEf)
str(spats_TDGxEf$predTrTot)

BLUEs_df <- as.data.frame(spats_TDGxEf$predTrTot)
str(BLUEs_df)
head(BLUEs_df)

# Renommer les colonnes
colnames(BLUEs_df) <- c("LIGNEE", "BLUEs_angle", "seBLUE")

# Vérifier le résultat
head(BLUEs_df)

#Enrégistrer les BLUES dans output
save(BLUEs_df, file = "output/pheno/BLUE_Angle_sordrought_2025.RData")

#Extraire la liste des génotypes
genotypes_BLUE <- unique(BLUEs_df$LIGNEE)
head(genotypes_BLUE)  # Vérifie les 6 premiers génotypes

#Enrégistrer les BLUES dans output
save(genotypes_BLUE, file = "output/geno/liste_Angle_genotypes.RData")



#=================================================================================================================================

library(SpATS)

#SPAT


data_root_Mean_angle$LIGNEE <- factor(data_root_Mean_angle$LIGNEE)
data_root_Mean_angle$Rep <- factor(data_root_Mean_angle$Rep)
data_root_Mean_angle$Block <- factor(data_root_Mean_angle$Block)

m <- SpATS(response = "Mean_angle_geno", genotype = "LIGNEE",
           genotype.as.random = TRUE,
           spatial = ~SAP(X, Y, nseg = c(20,20)),
           fixed = NULL, random = '~ Rep + Rep:Block',
           data = data_root_Mean_angle,
           control = list(maxit = 50, tolerance = 1e-06, monitoring = 1))



plot(m)

h2_SpATS <- getHeritability(m)

# BLUPs
pred <- predict(m, which = 'LIGNEE')
BLUP_SpATS <- pred$predicted.values
BLUP_SpATS <- pred[, c("LIGNEE", "predicted.values")]
colnames(BLUP_SpATS)[2] <- "BLUP_SpATS"

# compare lme4 and SpATS BLUP estimation
d_BLUP <- merge(BLUP_lme4, BLUP_SpATS, by = "LIGNEE")
plot(x = d_BLUP$BLUP_lme4, y = d_BLUP$BLUP_SpATS)

library(lme4)


#========================================================================================================================================================

THALLES



library(tidyr)
library(dplyr)
library(ggplot2)
library(agricolae)
library(multcompView)
library(statgenSTA)
library(dplyr)
library(readxl)


# setwd("C:/Users/2025an002/Desktop/Sordrought_2025/Aida_GWAS/Analyse_root_angle/data")
# list.files() 
data <- read_excel("data/pheno/data_angle_sordrought_2025.xlsx", sheet = 2)

geno_id <- unique(data$LIGNEE)
cross_id <- substr(x = data$LIGNEE, 1, 4)
table(cross_id)

# option add cross
# data$cross <- cross_id

str(data)

data$Thalles <- as.numeric(data$Thalles) 
#data <- na.omit(data)


#L? on obtient pour chaque g?notype, les donn?es ab?rentes pour l'ensemble des r?p?titions
data_Thalles <- data %>% 
  ggplot(aes(x = LIGNEE, y = Thalles)) +
  geom_boxplot(aes(x = LIGNEE, y = Thalles), outlier.shape = 8, outlier.color = "red", outlier.size = 3) +  
  geom_jitter(width = 0.2, alpha = 0.7, color = "blue") + 
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
data_Thalles


####### Fonction remove outliers
remove_outliers <- function(data, column) {
  Q1 <- quantile(data[[column]], 0.25, na.rm = TRUE)
  Q3 <- quantile(data[[column]], 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  
  data_filtered <- data[data[[column]] >= lower_bound & data[[column]] <= upper_bound, ]
  return(data_filtered)
}

# Loop through each Genotype and remove outliers
Genotypes <- unique(data$LIGNEE)
data_Thalles_cleaned <- data.frame()  # Empty dataframe to store results

for (gen in Genotypes) {
  subset_data <- data %>% filter(LIGNEE == gen)  # Subset for each Genotype
  cleaned_data <- remove_outliers(subset_data, "Thalles")  # Apply outlier removal
  data_Thalles_cleaned <- rbind(data_Thalles_cleaned, cleaned_data)  # Append cleaned data
}

print(data_Thalles_cleaned)

#write.table(data_Root_angle_clean, "data_Root_angle_clean.txt", sep = ";", row.names = FALSE, col.names = TRUE)
#Ici les outliers sont ?limin?s

data_Thalles_wo_Out <- data_Thalles_cleaned %>% 
  ggplot(aes(x = LIGNEE, y = Thalles)) +
  geom_boxplot(aes(x = LIGNEE, y = Thalles), outlier.shape = 8, outlier.color = "red", outlier.size = 3) +  
  geom_jitter(width = 0.2, alpha = 0.7, color = "blue") + 
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
data_Thalles_wo_Out


cat("Avant :", nrow(data), "\n")
cat("Après :", nrow(data_Thalles_cleaned), "\n")
cat("Supprimées :", nrow(data) - nrow(data_Thalles_cleaned), "\n")


#Avant : 2112 
#Après : 1988 
#Supprimées : 124 

##### Compute mean
# on calcule la moyenne des trois plants pour chaque g?notype. Chaque point repr?sente une r?p?tition

data_Thalles_wo_Out_means <- data_Thalles_cleaned %>%
  group_by(Y, X, X_Y, Block, Rep, Plot, LIGNEE) %>%
  summarise(Thalles = mean(Thalles), .groups = "drop")



data_Thalles_c <- data_Thalles_wo_Out_means %>% 
  ggplot(aes(x = LIGNEE, y = Thalles)) +
  geom_boxplot(aes(x = LIGNEE, y = Thalles), outlier.shape = 8, outlier.color = "red", outlier.size = 3) +  
  geom_jitter(width = 0.2, alpha = 0.7, color = "blue") + 
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
data_Thalles_c

########################################################
############# Correction StatGen#######################
################ On raw data###########################
########################################################

data_root_Mean_Thalles <- data_Thalles_wo_Out_means

data_root_Mean_Thalles$rowId <- as.factor(data_root_Mean_Thalles$Y)
data_root_Mean_Thalles$colId <- as.factor(data_root_Mean_Thalles$X)
data_root_Mean_Thalles$LIGNEE <- as.character(data_root_Mean_Thalles$LIGNEE)
data_root_Mean_Thalles$Thalles <- as.numeric(data_root_Mean_Thalles$Thalles)

data_root_Mean_Thalles <- data_root_Mean_Thalles %>%
  filter(
    is.finite(Thalles),
    !is.na(X),
    !is.na(Y)
  )

data_root_Mean_Thalles_R <- data_root_Mean_Thalles %>%
  dplyr::select(LIGNEE, X, Y, rowId, colId, Block, Rep, Plot, Thalles)



data_root_Mean_Thalles_R <- createTD(data = data_root_Mean_Thalles_R,
                                   genotype = "LIGNEE",
                                   rowCoord = "Y",
                                   colCoord = "X",
                                   rowId = "rowId",
                                   colId = "colId",
                                   subBlock = "Block",
                                   repId = "Rep",
                                   plotId = "Plot")

### Geno random
spaMod_data_Thalles_R <- fitTD(TD = data_root_Mean_Thalles_R,
                             design = "res.ibd",
                             traits = "Thalles",
                             what = "random")

summary(spaMod_data_Thalles_R)

## Create spatial plots of the results.
plot(spaMod_data_Thalles_R, plotType = "spatial", spaTrend = ("percentage"))

## Detect outliers in the standardized residuals of the fitted model.
spats_outliers <- outlierSTA(STA = spaMod_data_Thalles_R, traits = "Thalles", what = "random")
spats_outliers

## Extract all available statistics from the fitted model.
spats_extr <- extractSTA(spaMod_data_Thalles_R, what = "heritability")
spats_extr


#SUPPRESSION DES OUTLIERS
data_root_Mean_Thalles_R <- data_root_Mean_Thalles %>%
  dplyr::select(LIGNEE, X, Y, rowId, colId, Block, Rep, Plot, Thalles)

data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="459")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="640")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="557")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="484")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="397")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="485")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="318")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="57")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="63")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="593")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="696")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="521")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="1")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="284")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="25")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="492")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="145")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="221")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="406")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="121")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="487")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="486")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="161")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="191")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="573")
data_root_Mean_Thalles_R <- data_root_Mean_Thalles_R %>% filter(Plot!="348")


data_root_Mean_Thalles_R <- createTD(data = data_root_Mean_Thalles_R,
                                   genotype = "LIGNEE",
                                   rowCoord = "Y",
                                   colCoord = "X",
                                   rowId = "rowId",
                                   colId = "colId",
                                   subBlock = "Block",
                                   repId = "Rep",
                                   plotId = "Plot")

### Geno random
spaMod_data_Thalles_R <- fitTD(TD = data_root_Mean_Thalles_R,
                             design = "res.ibd",
                             traits = "Thalles",
                             what = "random")

summary(spaMod_data_Thalles_R)

## Create spatial plots of the results.
plot(spaMod_data_Thalles_R, plotType = "spatial", spaTrend = ("percentage"))

## Detect outliers in the standardized residuals of the fitted model.
spats_outliers <- outlierSTA(STA = spaMod_data_Thalles_R, traits = "Thalles", what = "random")
spats_outliers

## Extract all available statistics from the fitted model.
spats_extr <- extractSTA(spaMod_data_Thalles_R, what = "heritability")
spats_extr


### Geno fixed
data_root_Mean_Thalles_F <- data_root_Mean_Thalles %>%
  select(LIGNEE, X, Y, rowId, colId, Block, Rep, Plot, Thalles)


data_root_Mean_Thalles_F <- createTD(data = data_root_Mean_Thalles_F,
                                   genotype = "LIGNEE",
                                   rowCoord = "Y",
                                   colCoord = "X",
                                   rowId = "rowId",
                                   colId = "colId",
                                   subBlock = "Block",
                                   repId = "Rep",
                                   plotId = "Plot")

spaMod_Biom_F <- fitTD(TD = data_root_Mean_Thalles_F,
                       design = "res.ibd",
                       traits = "Thalles",
                       what = "fixed")

summary(spaMod_Biom_F)

## Create spatial plots of the results.
plot(spaMod_Biom_F, plotType = "spatial", spaTrend = ("percentage"))

## Detect outliers in the standardized residuals of the fitted model.
spats_outliers <- outlierSTA(STA = spaMod_Biom_F, traits = "Thalles", what = "fixed")
spats_outliers

#extr_lme4 <- extractSTA(spaMod, what=c("sMEANs"),  restoreColNames = TRUE, keep=c("repId","plotId","subBlock","rowCoord","colCoord"))
spats_TDGxEf <- STAtoTD(STA = spaMod_Biom_F, what = c("BLUEs", "seBLUEs"))
spats_TDGxEf




#Suppression des outliers

data_root_Mean_Thalles_F <- data_root_Mean_Thalles %>%
  select(LIGNEE, X, Y, rowId, colId, Block, Rep, Plot, Thalles)

data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="11")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="459")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="38")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="484")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="397")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="485")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="160")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="245")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="557")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="51")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="144")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="318")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="521")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="522")

data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="121")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="311")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="406")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="683")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="360")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="191")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="57")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="145")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="678")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="507")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="161")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="611")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="186")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="640")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="486")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="325")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="538")
data_root_Mean_Thalles_F <- data_root_Mean_Thalles_F %>% filter(Plot!="284")

data_root_Mean_Thalles_F <- createTD(data = data_root_Mean_Thalles_F,
                                   genotype = "LIGNEE",
                                   rowCoord = "Y",
                                   colCoord = "X",
                                   rowId = "rowId",
                                   colId = "colId",
                                   subBlock = "Block",
                                   repId = "Rep",
                                   plotId = "Plot")

spaMod_Biom_F <- fitTD(TD = data_root_Mean_Thalles_F,
                       design = "res.ibd",
                       traits = "Thalles",
                       what = "fixed")

summary(spaMod_Biom_F)

## Create spatial plots of the results.
plot(spaMod_Biom_F, plotType = "spatial", spaTrend = ("percentage"))

## Detect outliers in the standardized residuals of the fitted model.
spats_outliers <- outlierSTA(STA = spaMod_Biom_F, traits = "Thalles", what = "fixed")
spats_outliers
#extr_lme4 <- extractSTA(spaMod, what=c("sMEANs"),  restoreColNames = TRUE, keep=c("repId","plotId","subBlock","rowCoord","colCoord"))
spats_TDGxEf <- STAtoTD(STA = spaMod_Biom_F, what = c("BLUEs", "seBLUEs"))
spats_TDGxEf

#Write blups 
write.table(spats_TDGxEf, "BLUEs_Mean_Thalles_geno.txt", sep = ";", row.names = FALSE, col.names = TRUE)
write.csv(spats_TDGxEf, "BLUEs_Mean_Thalles_geno.csv")

write.csv(spats_TDGxEf, "data/pheno/BLUEs_Mean_Thalles_geno.csv")

#Pour enrégistrer les BLUES
str(spats_TDGxEf)
str(spats_TDGxEf$predTrTot)

BLUEs_df <- as.data.frame(spats_TDGxEf$predTrTot)
str(BLUEs_df)
head(BLUEs_df)

# Renommer les colonnes
colnames(BLUEs_df) <- c("LIGNEE", "BLUEs_Thalles", "seBLUE")

# Vérifier le résultat
head(BLUEs_df)

#Enrégistrer les BLUES dans output
save(BLUEs_df, file = "output/pheno/BLUE_Thalles_sordrought_2025.RData")


#Extraire la liste des génotypes
genotypes_BLUE <- unique(BLUEs_df$LIGNEE)
head(genotypes_BLUE)  # Vérifie les 6 premiers génotypes

#Enrégistrer les BLUES dans output
save(genotypes_BLUE, file = "output/geno/liste_Thalles_genotypes.RData")


#Pour les 2 fichier
# Charger les BLUEs pour Thalles
load("output/pheno/BLUE_Thalles_sordrought_2025.RData")
BLUEs_Thalles <- BLUEs_df

# Charger les BLUEs pour Angle
load("output/pheno/BLUE_Angle_sordrought_2025.RData")
BLUEs_Angle <- BLUEs_df

# Génotypes uniques pour Thalles
genotypes_Thalles <- unique(BLUEs_Thalles$LIGNEE)

# Génotypes uniques pour Angle
genotypes_Angle <- unique(BLUEs_Angle$LIGNEE)

# Liste complète des génotypes (sans doublons)
all_genotypes <- union(genotypes_Thalles, genotypes_Angle)

save(all_genotypes, file = "output/geno/liste_genotypes.RData")

