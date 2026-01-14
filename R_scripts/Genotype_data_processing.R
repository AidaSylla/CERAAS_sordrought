######################
# BCNAM genetic data #
######################
#Library----
library(readxl)
#Charger données----

load("data/par_per_cross/par_per_cross_GR.RData")
ppc_GR <- par_per_cross
# load("data/par_per_cross/par_per_cross_KK.RData")
# ppc_KK <- par_per_cross
rm(par_per_cross)

# les données de Aida provienne des croisements BC20, BC40, BC42

# source des données génétiques (dataverse)
# https://dataverse.cirad.fr/dataset.xhtml?persistentId=doi:10.18167/DVN1/TZVGLS
# mais il manque les données génétiques "brutes"

# ouvrir les données génétiques
load(file = 'data/geno/GR/geno_012_MAF.RData')
dim(geno_012)

geno_012_red <- geno_012[1:100, 1:100]

geno_012[1:10, 1:10]

# carte génétique
load(file = "data/map/Global_map.RData")

# identifier les lignées de la sous population: croisements BC20, BC40, BC42
data <- read_excel("data/pheno/data_angle_sordrought_2025.xlsx", sheet = 3)
X <- data$LIGNEE
geno_id <- unique(X)

#Subset des lignées "geno_id" présentes dans "geno_012" 
rownames(geno_012)
condition <- rownames(geno_012) %in% geno_id
geno <- geno_012[condition, ]

#Sauvegarde du fichier geno----
save(geno,file ="data/geno/GR/geno_subset_BC20_40_42.RData")







