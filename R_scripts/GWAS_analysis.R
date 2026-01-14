##### GWAS ANALYSIS #######
## ----chargement library GWAS-----------------------------------------------------------
library(tidyr)
library(dplyr)
library(rrBLUP)
library(AGHmatrix)
library(popkin)
library(qqman)

## ----chargement données GWAS---------------

load('./data/geno/GR/geno_subset_BC20_40_42.RData')
pheno <- read.csv('./data/pheno/BLUE/BLUEs_Mean_angle_geno.csv')
load('./data/map/Global_map.RData')

# subset des données génétique: 5K mks
# set.seed(5834)
# mk_sel <- sort(sample(x = 1:ncol(geno), size = 5000))
# geno <- geno[, mk_sel]
# map <- map[mk_sel, ]

#Transformation du fichier pheno ----
colnames(pheno) <- c("Plot", "id", "trait", "seBLUES")
pheno$cross <- sub("-.*", "", pheno$id)

#pheno2 <- pheno[, c(FALSE, TRUE, TRUE, FALSE, TRUE)]
pheno <- pheno[, c('id', 'trait',"cross")]
rownames(pheno) <- pheno$id

#Sélectionner les génotypes qui ont des données génétiques et phénotypiques
geno_com <- intersect(rownames(geno), pheno$id)
pheno <-  pheno[geno_com,]
geno <- geno[geno_com, ]

#---trait distribution----------------------------------------------------------------
hist(pheno$trait)

## -----repartition SNP-----------------------------------------------------------------

tapply(X = map$chr, INDEX = map$chr, FUN = length)

##--------Kinship matrix visualisation---------------------------------------------------

K_VR <- Gmatrix(SNPmatrix = geno, method = "VanRaden")

##--------Kinship matrix visualisation---------------------------------------------------
#plot_popkin(kinship = K_VR, labs = pheno$cross, width= 20)
cross_id <- substr(x = rownames(K_VR),start = 1,stop = 4)

plot_popkin(kinship = K_VR, labs = cross_id)


## formatage données génotypique
g_data <- data.frame(map[, -c(3)], t(geno), check.names = FALSE)

## ----classical GWAS--------------------------------------------------------------------
# Calcul du modèle
Q_GWAS <- rrBLUP::GWAS(pheno = pheno, fixed = 'cross', geno = g_data, K = K_VR,
                       plot = FALSE)

par(mfrow = c(1, 1))

# QQ-plot
p_val <- 10^(-Q_GWAS$trait)
qq(p_val)

# manhattan plot
d_man <- data.frame(CHR = map$chr, BP = map$phy_pos, P = p_val,
                    SNP = map$mk.names)
qqman::manhattan(x = d_man)

# selection QTLs

# seuil de sélection (conservateur)
bonfer_thre <- -log10(0.05/nrow(map))

QTL_bon <- Q_GWAS[Q_GWAS$trait > bonfer_thre, ]

QTL_log3 <- Q_GWAS[Q_GWAS$trait > 3, ]
QTL_log3 <- Q_GWAS[Q_GWAS$trait > 3.5, ]


#########################################################################################################################33
############################### Angle Racinaire #############################################################################333
###################################################################################################################
## ----chargement library GWAS-----------------------------------------------------------
library(tidyr)
library(dplyr)
library(rrBLUP)
library(AGHmatrix)
library(popkin)
library(qqman)

## ----chargement données GWAS---------------

load('geno/GR/geno_subset_BC20_40_42.RData')
pheno <- read.csv('pheno/BLUE/BLUEs_Mean_angle_geno.csv')
load('map/Global_map.RData')

# subset des données génétique: 5K mks
# set.seed(5834)
# mk_sel <- sort(sample(x = 1:ncol(geno), size = 5000))
# geno <- geno[, mk_sel]
# map <- map[mk_sel, ]

#Transformation du fichier pheno ----
colnames(pheno) <- c("Plot", "id", "root_angle", "seBLUES")
pheno$cross <- sub("-.*", "", pheno$id)

#pheno2 <- pheno[, c(FALSE, TRUE, TRUE, FALSE, TRUE)]
pheno <- pheno[, c('id', 'root_angle',"cross")]
rownames(pheno) <- pheno$id

#Sélectionner les génotypes qui ont des données génétiques et phénotypiques
geno_com <- intersect(rownames(geno), pheno$id)
pheno <-  pheno[geno_com,]
geno <- geno[geno_com, ]

#---root_angle distribution----------------------------------------------------------------
hist(pheno$root_angle)

## -----repartition SNP-----------------------------------------------------------------

tapply(X = map$chr, INDEX = map$chr, FUN = length)

##--------Kinship matrix visualisation---------------------------------------------------

K_VR <- Gmatrix(SNPmatrix = geno, method = "VanRaden")

##--------Kinship matrix visualisation---------------------------------------------------
#plot_popkin(kinship = K_VR, labs = pheno$cross, width= 20)
cross_id <- substr(x = rownames(K_VR),start = 1,stop = 4)

plot_popkin(kinship = K_VR, labs = cross_id)


## formatage données génotypique
g_data <- data.frame(map[, -c(3)], t(geno), check.names = FALSE)

## ----classical GWAS--------------------------------------------------------------------
# Calcul du modèle
Q_GWAS <- rrBLUP::GWAS(pheno = pheno, fixed = 'cross', geno = g_data, K = K_VR,
                       plot = FALSE)

par(mfrow = c(1, 1))

# QQ-plot
p_val <- 10^(-Q_GWAS$root_angle)
qq(p_val)

# manhattan plot
d_man <- data.frame(CHR = map$chr, BP = map$phy_pos, P = p_val,
                    SNP = map$mk.names)
qqman::manhattan(x = d_man)

# selection QTLs

# seuil de sélection (conservateur)
bonfer_thre <- -log10(0.05/nrow(map))

QTL_bon <- Q_GWAS[Q_GWAS$root_angle > bonfer_thre, ]

QTL_log3 <- Q_GWAS[Q_GWAS$root_angle > 3, ]
QTL_log3 <- Q_GWAS[Q_GWAS$root_angle > 3.5, ]




#######################################################################################################################################################################################################
######################### THALLES ############################################################################################
#######################################################################################################################################################################################################

## ----chargement library GWAS-----------------------------------------------------------
library(tidyr)
library(dplyr)
library(rrBLUP)
library(AGHmatrix)
library(popkin)
library(qqman)

## ----chargement données GWAS---------------

load('geno/GR/geno_subset_BC20_40_42.RData')
pheno <- read.csv('pheno/BLUE/BLUEs_Mean_Thalles_geno.csv')
load('map/Global_map.RData')

# subset des données génétique: 5K mks
# set.seed(5834)
# mk_sel <- sort(sample(x = 1:ncol(geno), size = 5000))
# geno <- geno[, mk_sel]
# map <- map[mk_sel, ]

#Transformation du fichier pheno ----
colnames(pheno) <- c("Plot", "id", "Thalles", "seBLUES")
pheno$cross <- sub("-.*", "", pheno$id)

#pheno2 <- pheno[, c(FALSE, TRUE, TRUE, FALSE, TRUE)]
pheno <- pheno[, c('id', 'Thalles',"cross")]
rownames(pheno) <- pheno$id

#Sélectionner les génotypes qui ont des données génétiques et phénotypiques
geno_com <- intersect(rownames(geno), pheno$id)
pheno <-  pheno[geno_com,]
geno <- geno[geno_com, ]

#---Thalles distribution----------------------------------------------------------------
hist(pheno$Thalles)

## -----repartition SNP-----------------------------------------------------------------

tapply(X = map$chr, INDEX = map$chr, FUN = length)

##--------Kinship matrix visualisation---------------------------------------------------

K_VR <- Gmatrix(SNPmatrix = geno, method = "VanRaden")

##--------Kinship matrix visualisation---------------------------------------------------
#plot_popkin(kinship = K_VR, labs = pheno$cross, width= 20)
cross_id <- substr(x = rownames(K_VR),start = 1,stop = 4)

plot_popkin(kinship = K_VR, labs = cross_id)


## formatage données génotypique
g_data <- data.frame(map[, -c(3)], t(geno), check.names = FALSE)

## ----classical GWAS--------------------------------------------------------------------
# Calcul du modèle
Q_GWAS <- rrBLUP::GWAS(pheno = pheno, fixed = 'cross', geno = g_data, K = K_VR,
                       plot = FALSE)

par(mfrow = c(1, 1))

# QQ-plot
p_val <- 10^(-Q_GWAS$Thalles)
qq(p_val)

# manhattan plot
d_man <- data.frame(CHR = map$chr, BP = map$phy_pos, P = p_val,
                    SNP = map$mk.names)
qqman::manhattan(x = d_man)

# selection QTLs

# seuil de sélection (conservateur)
bonfer_thre <- -log10(0.05/nrow(map))

QTL_bon <- Q_GWAS[Q_GWAS$Thalles > bonfer_thre, ]

QTL_log3 <- Q_GWAS[Q_GWAS$Thalles > 3, ]
QTL_log3 <- Q_GWAS[Q_GWAS$Thalles > 4, ]



