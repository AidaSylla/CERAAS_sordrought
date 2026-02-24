###############################
# mppData object construction #
###############################

library(mppR)

# determine from which recurrent parent the cross come from ----

# list of genotpes
load(file = "output/geno/liste_Angle_genotypes.RData")

geno_id_angle <- genotypes_BLUE

load(file = "output/geno/liste_Thalles_genotypes.RData")

geno_id_thalle <- genotypes_BLUE

load(file = "output/geno/liste_genotypes.RData")

geno_id <- unique(union(union(as.character(geno_id_angle),
                        as.character(geno_id_thalle)), as.character(all_genotypes)))

cross_ind <- substr(x = geno_id, 1, 4)

# par_per_cross
load(file = "data/par_per_cross/par_per_cross_KK.RData")
ppc_KK <- par_per_cross

load(file = "data/par_per_cross/par_per_cross_GR.RData")
ppc_GR <- par_per_cross

par_per_cross <- rbind(ppc_KK, ppc_GR)

ppc_data <- par_per_cross[par_per_cross$cross %in% cross_ind, ] 

unique(ppc_data$Par1)

# Only data from the Grinkan sub-population.

# Subset from the Grinkan global mppData ----

load(file = "data/mppData/raw_data/GR/mppData.RData")

# subset: keep only the genotype present
sum(geno_id %in% mppData$geno.id)
geno_com <- intersect(geno_id, mppData$geno.id)

mppData <- subset(mppData, gen.list = geno_com)
mppData$geno.off <- NULL

# replace the phenotype data ----

# load data
load(file = "output/pheno/BLUE_Angle_sordrought_2025.RData")
BLUE_angle <- BLUEs_df
rownames(BLUE_angle) <- BLUE_angle$LIGNEE

load(file = "output/pheno/BLUE_Thalles_sordrought_2025.RData")
BLUE_thalle <- BLUEs_df
rownames(BLUE_thalle) <- BLUE_thalle$LIGNEE

ref_geno <- data.frame(geno = mppData$geno.id)
ref_geno$angle <- BLUE_angle[ref_geno$geno, ]$BLUEs_angle
ref_geno$thalle <- BLUE_thalle[ref_geno$geno, ]$BLUEs_Thalles

identical(ref_geno$geno, mppData$geno.id)

pheno <- ref_geno
rownames(pheno) <- pheno$geno
pheno <- as.matrix(pheno[, -1])

mppData$pheno <- pheno

# save a copy of the subseted mppData object ----
#save(mppData, file = "data/mppData/mppData.RData")
#Vincent m'a demandé de ne pas exécuter la ligne 75


#========================================================================================================
#========================================================================================================
# QTL detection (test) ----
load(file = "data/mppData/mppData.RData")

#perm <- mpp_perm(mppData, trait = "angle", Q.eff = "par")
#perm$threshold
perm <- mpp_perm(
  mppData = mppData,
  trait   = "angle",
  Q.eff   = "par",
  N       = 1000
)

perm$threshold

save(perm, file = "perm_angle_1000.RData")

#Simple interval mapping (SIM)
#SIM <- mpp_SIM(
#  mppData = mppData,
#  trait   = "angle",
#  Q.eff   = "par",
#  threshold = perm$threshold,
#  plot.gen.eff = TRUE
#)

#Accélérer avec plusieurs cœurs 
perm <- mpp_perm(
  mppData = mppData,
  trait   = "angle",
  Q.eff   = "par",
  N       = 1000,
  n.cores = 4
)

perm$threshold

#Pour sauver le resultat
save(perm, file = "results/perm_1000_angle_Qeff_par.RData")

load("results/perm_1000_angle_Qeff_par.RData")

#thr_95 <- perm$quantiles["95%"]
#thr_95
# autre méthode
thr_95 <- perm$quantiles["95%"]
saveRDS(thr_95, file = "results/seuil_perm_angle_95.rds")


#thr_95 <- readRDS("results/seuil_perm_angle_95.rds")
#thr_95

#Extraction correcte du seuil 95 %
thr_95 <- perm$q.val
thr_95

thr_95 <- unname(perm$q.val)
thr_95

saveRDS(thr_95, file = "results/seuil_perm_1000_angle_95.rds")

thr_95 <- readRDS("results/seuil_perm_1000_angle_95.rds")
thr_95


save(perm, file = "results/perm_1000_angle_Qeff_par.RData")

load("results/perm_1000_angle_Qeff_par.RData")
thr_95 <- unname(perm$q.val)

#SIMPLE INTERVAL MAPPING (SIM)
library(mppR)
SIM <- mpp_SIM(mppData = mppData, trait = "angle", Q.eff = "par",
               plot.gen.eff = TRUE)
plot(SIM)
plot(SIM, gen.eff = TRUE, Q.eff = "par", mppData = mppData)

# SIM seuil Permitation 1000
QTL <- QTL_select(
  Qprof = SIM,
  threshold = thr_95
)
plot(SIM, threshold = thr_95)
plot(SIM, threshold = thr_95, gen.eff = TRUE, Q.eff = "par", mppData = mppData)

# Vérifier le QTL sélectionné
print(QTL)

# A compléter et étendre dans un nouveau script

#Cofactors selection
cofactors <- QTL_select(Qprof = SIM, threshold = thr_95, window = 50)

print(cofactors)

#Composite interval mapping (CIM)
### CIM (sans **threshold** dans l’appel)
CIM <- mpp_CIM(
  mppData = mppData,
  trait = "angle",
  cofactors = cofactors,
  window = 20,
  plot.gen.eff = TRUE,
  n.cores = 4
)

### sélection des QTL à partir du profil CIM
QTL <- QTL_select(
  Qprof = CIM,
  threshold = thr_95,
  window = 20
)
plot(CIM, threshold = thr_95)

str(QTL)

# Vérifier le QTL sélectionné
print(QTL)

# Estimer les effets génétiques pour tous les QTL détectés
effects <- QTL_gen_effects(
  mppData = mppData,
  trait   = "angle",
  QTL     = QTL,     # liste de QTL
  Q.eff   = "par",   # parental effect
  sum_zero = TRUE    # optionnelle
)

# Vérifier les résultats
effects$Qeff  # liste avec effets par QTL
effects$tab.Qeff  # tableau récapitulatif avec chromosome, position et effets

#QTLs R2 calculation
QR2 <- QTL_R2(mppData = mppData, trait = "angle", QTL = QTL)

#la contribution des QTL
QR2$glb.adj.R2

# la contribution de chaque QTL
QR2$part.adj.R2.diff


#=================================

QTL_eff <- QTL_gen_effects(
  mppData = mppData,
  trait   = "angle",
  QTL     = QTL,
  Q.eff   = "par"
)
#Afficher les effets
summary(QTL_eff)

plot(QTL_eff)

str(CIM, max.level = 2)






# Exemple pour le QTL sur chr 1 à 17.7193 cM
qtl_chr <- QTL$chr
qtl_pos <- QTL$pos.cM

# Subset du chromosome
chr_prof <- subset(CIM, chr == qtl_chr)

# Trouver LOD max proche du QTL
LODmax <- max(chr_prof$log10pval[abs(chr_prof$pos.cM - qtl_pos) < 1])  # ±1 cM autour du QTL

# Définir le seuil pour IC 1-LOD
LOD_IC <- LODmax - 1

# Positions où LOD >= LOD_IC
positions_IC <- chr_prof$pos.cM[chr_prof$log10pval >= LOD_IC]

# Intervalle de confiance
IC_left  <- min(positions_IC)
IC_right <- max(positions_IC)

IC <- data.frame(chr = qtl_chr, QTL_pos = qtl_pos, left = IC_left, right = IC_right)
IC




















#========================================================================================================
#========================================================================================================
##============================THALLES====================================================================
#========================================================================================================
#========================================================================================================
# QTL detection (test) ----
load(file = "data/mppData/mppData.RData")

library(mppR)
#perm <- mpp_perm(mppData, trait = "angle", Q.eff = "par")
perm <- mpp_perm(
  mppData = mppData,
  trait   = "thalle",   # en minuscules
  Q.eff   = "par",
  N       = 1000
)

perm$threshold

save(perm, file = "perm_angle_1000.RData")


#Accélérer avec plusieurs cœurs 
perm <- mpp_perm(
  mppData = mppData,
  trait   = "thalle",
  Q.eff   = "par",
  N       = 1000,
  n.cores = 4
)

perm$threshold

#Pour sauver le resultat
save(perm, file = "results/perm_1000_angle_Qeff_par.RData")

load("results/perm_1000_angle_Qeff_par.RData")

#thr_95 <- perm$quantiles["95%"]
#thr_95
# autre méthode
thr_thalle_95 <- perm$quantiles["95%"]
saveRDS(thr_thalle_95, file = "results/seuil_perm_angle_95.rds")


#thr_95 <- readRDS("results/seuil_perm_angle_95.rds")
#thr_95

#Extraction correcte du seuil 95 %
thr_thalle_95 <- perm$q.val
thr_thalle_95

thr_thalle_95 <- unname(perm$q.val)
thr_thalle_95

saveRDS(thr_thalle_95, file = "results/seuil_perm_1000_angle_95.rds")

thr_thalle_95 <- readRDS("results/seuil_perm_1000_angle_95.rds")
thr_thalle_95


save(perm, file = "results/perm_1000_thalle_Qeff_par.RData")

load("results/perm_1000_thalle_Qeff_par.RData")
thr_thalle_95 <- unname(perm$q.val)

#SIMPLE INTERVAL MAPPING (SIM)
library(mppR)
SIM <- mpp_SIM(mppData = mppData, trait = "thalle", Q.eff = "par",
               plot.gen.eff = TRUE)
plot(SIM)
plot(SIM, gen.eff = TRUE, Q.eff = "par", mppData = mppData)

# SIM seuil Permitation 1000
QTL <- QTL_select(
  Qprof = SIM,
  threshold = thr_thalle_95
)
plot(SIM, threshold = thr_thalle_95)
plot(SIM, threshold = thr_thalle_95, gen.eff = TRUE, Q.eff = "par", mppData = mppData)

# Vérifier le QTL sélectionné
print(QTL)

# A compléter et étendre dans un nouveau script

#Cofactors selection
cofactors <- QTL_select(Qprof = SIM, threshold = thr_thalle_95, window = 50)

print(cofactors)

#Composite interval mapping (CIM)
### CIM (sans **threshold** dans l’appel)
CIM <- mpp_CIM(
  mppData = mppData,
  trait = "thalle",
  cofactors = cofactors,
  window = 20,
  plot.gen.eff = TRUE,
  n.cores = 4
)

### sélection des QTL à partir du profil CIM
QTL <- QTL_select(
  Qprof = CIM,
  threshold = thr_thalle_95,
  window = 20
)
plot(CIM, threshold = thr_thalle_95)

str(QTL)

# Vérifier le QTL sélectionné
print(QTL)

# Estimer les effets génétiques pour tous les QTL détectés
effects <- QTL_gen_effects(
  mppData = mppData,
  trait   = "angle",
  QTL     = QTL,     # liste de QTL
  Q.eff   = "par",   # parental effect
  sum_zero = TRUE    # optionnelle
)

# Vérifier les résultats
effects$Qeff  # liste avec effets par QTL
effects$tab.Qeff  # tableau récapitulatif avec chromosome, position et effets

#QTLs R2 calculation
QR2 <- QTL_R2(mppData = mppData, trait = "angle", QTL = QTL)

#la contribution des QTL
QR2$glb.adj.R2

# la contribution de chaque QTL
QR2$part.adj.R2.diff


#=================================

QTL_eff <- QTL_gen_effects(
  mppData = mppData,
  trait   = "angle",
  QTL     = QTL,
  Q.eff   = "par"
)
#Afficher les effets
summary(QTL_eff)

plot(QTL_eff)

str(CIM, max.level = 2)






# Exemple pour ton QTL sur chr 1 à 17.7193 cM
qtl_chr <- QTL$chr
qtl_pos <- QTL$pos.cM

# Subset du chromosome
chr_prof <- subset(CIM, chr == qtl_chr)

# Trouver LOD max proche du QTL
LODmax <- max(chr_prof$log10pval[abs(chr_prof$pos.cM - qtl_pos) < 1])  # ±1 cM autour du QTL

# Définir le seuil pour IC 1-LOD
LOD_IC <- LODmax - 1

# Positions où LOD >= LOD_IC
positions_IC <- chr_prof$pos.cM[chr_prof$log10pval >= LOD_IC]

# Intervalle de confiance
IC_left  <- min(positions_IC)
IC_right <- max(positions_IC)

IC <- data.frame(chr = qtl_chr, QTL_pos = qtl_pos, left = IC_left, right = IC_right)
IC