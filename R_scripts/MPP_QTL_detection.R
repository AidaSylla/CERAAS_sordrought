#####################
# MPP QTL detection #
#####################

# Install mppR (github)
# library(devtools)
# devtools::install_github(repo = "vincentgarin/mppR", ref = "master")

library(mppR)
#========================================================================================================
#========================================================================================================
# QTL detection (test) ----
load(file = "data/mppData/mppData.RData")


#===============================================
#CONTROLE QUALITE SNPs
#================================================
# Matrice SNP
geno <- mppData$geno.IBS

# 1. Missing rate par SNP
missing_snp <- colMeans(is.na(geno))

summary(missing_snp)

# garder SNP avec <=10% missing
keep_missing <- missing_snp <= 0.10

# 2. Calcul MAF pour codage 0/1/2
maf <- apply(geno, 2, function(x){
  
  p <- mean(x, na.rm = TRUE) / 2
  
  maf <- min(p, 1 - p)
  
  return(maf)
  
})

summary(maf)

# garder SNP MAF >= 5%
keep_maf <- maf >= 0.05

# 3. Filtre final
keep_SNP <- keep_missing & keep_maf


table(keep_SNP)

cat("Nombre SNP avant QC :", ncol(geno), "\n")
cat("Nombre SNP après QC :", sum(keep_SNP), "\n")


mppData <- mppData

mppData$geno.IBS <- mppData$geno.IBS[, keep_SNP]

dim(mppData$geno.IBS)


save(
  mppData,
  file = "data/mppData/mppData.RData"
)

#Après celà les Snps sont filtrés


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


load("threshold_angle_perm.RData")


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