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
setwd("C:/Users/2025an002/Desktop/Sordrought_2025/Aida_GWAS")
data <- read_excel("data/pheno/data_angle_sordrought_2025.xlsx", sheet = 3)

geno_id <- unique(data$LIGNEE)
cross_id <- substr(x = data$LIGNEE, 1, 4)
table(cross_id)

# option add cross
# data$cross <- cross_id



#ANALYSE SPATIALE AVEC STATGEN ET IQR
colnames(data)

data$Root_angle <- as.numeric(data$Root_angle) 

data1 <- data %>% 
  ggplot(aes(x = LIGNEE, y = Root_angle)) +
  geom_boxplot(aes(x = LIGNEE, y = Root_angle), outlier.shape = 8, outlier.color = "red", outlier.size = 3) +  
  geom_jitter(width = 0.2, alpha = 0.7, color = "blue") + 
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
data1



#Pour que les outliers ne soient pas montrés

data1 <- data %>% 
  ggplot(aes(x = LIGNEE, y = Root_angle)) +
  geom_boxplot(outlier.shape = NA) +  # ne montre pas les outliers
  geom_jitter(width = 0.2, alpha = 0.7, color = "blue") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

data1

####### Fonction remove outliers
####### Fonction remove outliers
#remove_outliers <- function(data, column) {
#  Q1 <- quantile(data[[column]], 0.25, na.rm = TRUE)
#  Q3 <- quantile(data[[column]], 0.75, na.rm = TRUE)
#  IQR <- Q3 - Q1
  
#  lower_bound <- Q1 - 1.5 * IQR
#  upper_bound <- Q3 + 1.5 * IQR
  
#  data_filtered <- data[data[[column]] >= lower_bound & data[[column]] <= upper_bound, ]
#  return(data_filtered)
# }

# Loop through each LIGNEE and remove outliers
#LIGNEEs <- unique(data$LIGNEE)
#data_angle_cleaned <- data.frame()  # Empty dataframe to store results

#for (gen in LIGNEEs) {
#  subset_data <- data %>% filter(LIGNEE == gen)  # Subset for each LIGNEE
#  cleaned_data <- remove_outliers(subset_data, "Root_angle")  # Apply outlier removal
#  data_angle_cleaned <- rbind(data_angle_cleaned, cleaned_data)  # Append cleaned data
# }

#print(data_angle_cleaned)

##### Compute mean
# on calcule la moyenne des trois plants pour chaque g?notype. Chaque point repr?sente une r?p?tition
data_angle_means <- data %>% 
  group_by(LIGNEE, Plot) %>%
  summarise(Root_angle = mean(Root_angle, na.rm = TRUE))

data_angle_c <- data_angle_means %>% 
  ggplot(aes(x = LIGNEE, y = Root_angle)) +
  geom_boxplot(aes(x = LIGNEE, y = Root_angle), outlier.shape = 8, outlier.color = "red", outlier.size = 3) +  
  geom_jitter(width = 0.2, alpha = 0.7, color = "blue") + 
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
data_angle_c


#==========================================================================================================================================================================
#=============== Correction StatGen =======================================================================================================================================
#================== On raw data ==============================================================================================================================================
#===========================================================================================================================================================================
setwd("C:/Users/2025an002/Desktop/Sordrought_2025/Aida_GWAS")
data_root_Mean_angle <- read_excel("data/pheno/data_angle_sorDrought_2025.xlsx", sheet = 4)

data_root_Mean_angle$rowId <- as.factor(data_root_Mean_angle$Y)
data_root_Mean_angle$colId <- as.factor(data_root_Mean_angle$X)
data_root_Mean_angle$LIGNEE <- as.character(data_root_Mean_angle$LIGNEE)
data_root_Mean_angle$Mean_angle_geno <- as.numeric(data_root_Mean_angle$Mean_angle_geno)



data_root_Mean_angle_R <- data_root_Mean_angle %>%
  select(LIGNEE, X, Y, rowId, colId, Block, Rep, Plot, Mean_angle_geno)

data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="562")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="412")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="500")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="597")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="352")


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
                             traits = "Mean_angle_geno",
                             what = "random")

summary(spaMod_data_angle_R)

## Create spatial plots of the results.
plot(spaMod_data_angle_R, plotType = "spatial", spaTrend = ("percentage"))

## Detect outliers in the standardized residuals of the fitted model.
spats_outliers <- outlierSTA(STA = spaMod_data_angle_R, traits = "Mean_angle_geno", what = "random")

## Extract all available statistics from the fitted model.
spats_extr <- extractSTA(spaMod_data_angle_R, what = "heritability")
spats_extr

# BLUES <- extract()

### Geno fixed
data_root_Mean_angle_F <- data_root_Mean_angle %>%
  select(LIGNEE, X, Y, rowId, colId, Block, Rep, Plot, Mean_angle_geno)

data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="104")
data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="690")
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
                       traits = "Mean_angle_geno",
                       what = "fixed")

summary(spaMod_Biom_F)

## Create spatial plots of the results.
plot(spaMod_Biom_F, plotType = "spatial", spaTrend = ("percentage"))

## Detect outliers in the standardized residuals of the fitted model.
spats_outliers <- outlierSTA(STA = spaMod_Biom_F, traits = "Mean_angle_geno", what = "fixed")

#extr_lme4 <- extractSTA(spaMod, what=c("sMEANs"),  restoreColNames = TRUE, keep=c("repId","plotId","subBlock","rowCoord","colCoord"))
spats_TDGxEf <- STAtoTD(STA = spaMod_Biom_F, what = c("BLUEs", "seBLUEs"))
spats_TDGxEf

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
save(BLUEs_df, file = "output/BLUE_Angle_sordrought_2025.RData")

#Extraire la liste des génotypes
genotypes_BLUE <- unique(BLUEs_df$LIGNEE)
head(genotypes_BLUE)  # Vérifie les 6 premiers génotypes

#Enrégistrer les BLUES dans output
save(genotypes_BLUE, file = "output/liste_Angle_genotypes.RData")



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


setwd("C:/Users/2025an002/Desktop/Sordrought_2025/CERAAS_sordrought/data")
# list.files() 
data <- read_excel("pheno/data_angle_sordrought_2025.xlsx", sheet = 3)

geno_id <- unique(data$LIGNEE)
cross_id <- substr(x = data$LIGNEE, 1, 4)
table(cross_id)

#ANALYSE SPATIALE AVEC STATGEN ET IQR
colnames(data)

data$Thalles <- as.numeric(data$Thalles) 


data1 <- data %>% 
  ggplot(aes(x = LIGNEE, y = Thalles)) +
  geom_boxplot(aes(x = LIGNEE, y = Thalles), outlier.shape = 8, outlier.color = "red", outlier.size = 3) +  
  geom_jitter(width = 0.2, alpha = 0.7, color = "blue") + 
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
data1



data1 <- data %>% 
  ggplot(aes(x = LIGNEE, y = Thalles)) +
  geom_boxplot(outlier.shape = NA) +  # ne montre pas les outliers
  geom_jitter(width = 0.2, alpha = 0.7, color = "blue") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

data1

#data_angle_mean=data_angle_means
##### Compute mean
# on calcule la moyenne des trois plants pour chaque g?notype. Chaque point repr?sente une r?p?tition
#J'ai mis na.rm=TRUE pour éviter qu'il supprime toutes les valeurs d'une parcelle s"il y'a une donnée manquante
data_angle_means <- data %>%
  group_by(LIGNEE, Plot) %>%
  summarise(Thalles = mean(Thalles, na.rm = TRUE))

data_angle_c <- data_angle_means %>% 
  ggplot(aes(x = LIGNEE, y = Thalles)) +
  geom_boxplot(aes(x = LIGNEE, y = Thalles), outlier.shape = 8, outlier.color = "red", outlier.size = 3) +  
  geom_jitter(width = 0.2, alpha = 0.7, color = "blue") + 
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
data_angle_c


#=====================================================================================================================
############# Correction StatGen#######################
################ On raw data###########################
########################################################
setwd("C:/Users/2025an002/Desktop/Sordrought_2025/CERAAS_sordrought")
data_root_Mean_angle <- read_excel("data/pheno/data_angle_sordrought_2025.xlsx", sheet = 4)

data_root_Mean_angle$rowId <- as.factor(data_root_Mean_angle$Y)
data_root_Mean_angle$colId <- as.factor(data_root_Mean_angle$X)
data_root_Mean_angle$LIGNEE <- as.character(data_root_Mean_angle$LIGNEE)
data_root_Mean_angle$Mean_Thalles_geno <- as.numeric(data_root_Mean_angle$Mean_Thalles_geno)



data_root_Mean_angle_R <- data_root_Mean_angle %>%
  select(LIGNEE, X, Y, rowId, colId, Block, Rep, Plot, Mean_Thalles_geno)

data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="191")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="371")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="459")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="639")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="557")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="484")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="485")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="318")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="406")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="63")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="590")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="593")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="696")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="521")

data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="640")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="25")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="396")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="492")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="494")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="57")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="692")

data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="573")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="404")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="322")

data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="182")
data_root_Mean_angle_R <- data_root_Mean_angle_R %>% filter(Plot!="145")

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
                             traits = "Mean_Thalles_geno",
                             what = "random")

summary(spaMod_data_angle_R)

## Create spatial plots of the results.
plot(spaMod_data_angle_R, plotType = "spatial", spaTrend = ("percentage"))

## Detect outliers in the standardized residuals of the fitted model.
spats_outliers <- outlierSTA(STA = spaMod_data_angle_R, traits = "Mean_Thalles_geno", what = "random")

## Extract all available statistics from the fitted model.
spats_extr <- extractSTA(spaMod_data_angle_R, what = "heritability")
spats_extr

### Geno fixed
data_root_Mean_angle_F <- data_root_Mean_angle %>%
  select(LIGNEE, X, Y, rowId, colId, Block, Rep, Plot, Mean_Thalles_geno)

data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="11")
data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="484")
data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="191")
data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="507")
data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="38")
data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="459")
data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="144")
data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="557")
data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="311")
data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="406")
data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="245")
data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="485")
data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="65")
data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="494")

data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="105")
data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="639")
data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="122")
data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="696")
data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="51")
data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="521")

data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="182")
data_root_Mean_angle_F <- data_root_Mean_angle_F %>% filter(Plot!="641")


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
                       traits = "Mean_Thalles_geno",
                       what = "fixed")

summary(spaMod_Biom_F)

## Create spatial plots of the results.
plot(spaMod_Biom_F, plotType = "spatial", spaTrend = ("percentage"))

## Detect outliers in the standardized residuals of the fitted model.
spats_outliers <- outlierSTA(STA = spaMod_Biom_F, traits = "Mean_Thalles_geno", what = "fixed")

#extr_lme4 <- extractSTA(spaMod, what=c("sMEANs"),  restoreColNames = TRUE, keep=c("repId","plotId","subBlock","rowCoord","colCoord"))
spats_TDGxEf <- STAtoTD(STA = spaMod_Biom_F, what = c("BLUEs", "seBLUEs"))
spats_TDGxEf

#Write blues 
#write.table(spats_TDGxEf, "BLUEs_Mean_Thalles_geno.txt", sep = ";", row.names = FALSE, col.names = TRUE)
#write.csv(spats_TDGxEf, "BLUEs_Mean_Thalles_geno.csv")

#write.csv(spats_TDGxEf, "pheno/BLUE/BLUEs_Mean_Thalles_geno.csv")


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

