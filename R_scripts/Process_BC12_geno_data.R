#============================================================================================================================
#La première partie est le script que Vincent m'a envoyé
#============================================================================================================================

######################
# BC12 genotype data #
######################

library(qtl)
library(dplyr)

# modify the map data
# load(file = "data/map/Global_map.RData")
# 
# geno_raw <- read.table(file = "data/geno/GR/BC12_final_RQTL_input.txt",
#                        fill = NA, sep = ",")
# 
# geno_raw[3, 4:ncol(geno_raw)] <- map$pos.cM
# 
# write.table(geno_raw, file = "data/geno/GR/BC12_final_RQTL_input.txt",
#             col.names = FALSE, row.names = FALSE, quote = FALSE, sep = ",")

cr_obj <- read.cross(format="csv", file = "data/geno/GR/BC12_final_RQTL_input.txt",
                    na.strings="NA", error.prob=0.05,
                    map.function="Kosambi", BC.gen=1, F.gen=3)

geno.image(cr_obj)

# modify the phenotype of the cross object to add your trait data
pheno_cr_obj <- cr_obj$pheno
geno_id_org <- pheno_cr_obj$id
pheno_cr_obj$id <- toupper(pheno_cr_obj$id)

# merge your phenotype
my_pheno <- data.frame(id = pheno_cr_obj$id,
                      angle = rnorm(n = nrow(pheno_cr_obj)))

pheno_cr_obj <- left_join(pheno_cr_obj, my_pheno)
pheno_cr_obj$id <- tolower(pheno_cr_obj$id)
pheno_cr_obj <- pheno_cr_obj[, -c(2:3)]

cr_obj$pheno <- pheno_cr_obj

SIM <- scanone(cross = cr_obj, pheno.col = 2, method = "hk")
plot(SIM, col = "blue")



# Afficher les phénotypes actuels
head(cr_obj$pheno)

# Afficher les noms des colonnes phénotypiques
colnames(cr_obj$pheno)



##################################################################################################################





























#============================================================================================================================
#EN BAS, C'est le script que j'ai fais pour intégréer le phénotype étudié (TR_RATE).
#============================================================================================================================


library(qtl)
library(dplyr)

---
  ### 🔹 Étape 1 : Charger l'objet cross (cr_obj)
  print("🔹 Chargement de l'objet cross...")
cr_obj <- read.cross(
  format = "csv",
  file = "data/geno/GR/BC12_final_RQTL_input.txt",
  na.strings = "NA",
  error.prob = 0.05,
  map.function = "Kosambi",
  BC.gen = 1,
  F.gen = 3
)

# Visualiser les génotypes (optionnel)
geno.image(cr_obj)

  ### 🔹 Étape 2 : Charger et nettoyer les données phénotypiques (TR_RATE)
  print("\n🔹 Chargement et nettoyage de TR_RATE...")

# Charger le fichier BLUEs_TR_RATE_2026.csv
pheno_data <- read.csv2("data/pheno/BLUEs_Pente.csv")

# Renommer la colonne "Lines" en "id"
colnames(pheno_data)[colnames(pheno_data) == "GENOTYPE"] <- "id"

# Convertir TR_RATE en numérique
pheno_data$BLUE_Slope <- as.numeric(pheno_data$BLUE_Slope)

# Nettoyer les identifiants : convertir en minuscules et remplacer "BC12-" par "bc12-"
pheno_data$id <- tolower(sub("^BC12-", "bc12-", pheno_data$id))

---
  ### 🔹 Étape 3 : Filtrer pheno_data pour ne garder que les individus présents dans cr_obj
  print("\n🔹 Filtrage de pheno_data...")

# Convertir les id de cr_obj$pheno en minuscules
cr_obj$pheno$id <- tolower(cr_obj$pheno$id)

# Garder seulement les individus de pheno_data qui sont dans cr_obj$pheno
pheno_data <- pheno_data[pheno_data$id %in% cr_obj$pheno$id, ]

# Vérifier le nombre d'individus après filtrage
print(paste("✅ Nombre d'individus dans pheno_data après filtrage :", nrow(pheno_data)))

---
  ### 🔹 Étape 4 : Fusionner TR_RATE avec cr_obj$pheno
  print("\n🔹 Fusion de TR_RATE avec cr_obj$pheno...")

# Fusionner les données
cr_obj$pheno <- left_join(cr_obj$pheno, pheno_data, by = "id")

# Vérifier que TR_RATE est bien ajouté
print("✅ Phénotypes après fusion :")
head(cr_obj$pheno)
print(paste("📊 Nombre de NA dans TR_RATE :", sum(is.na(cr_obj$pheno$TR_RATE))))

---
  ### 🔹 Étape 5 : Supprimer pheno_1 et pheno_2
  print("\n🔹 Suppression de pheno_1 et pheno_2...")

# Supprimer pheno_1 et pheno_2
cr_obj$pheno <- cr_obj$pheno %>% select(id, BLUE_Slope)

# Vérifier le résultat
print("✅ Structure de cr_obj$pheno après suppression :")
head(cr_obj$pheno)
colnames(cr_obj$pheno)

---
  ### 🔹 Étape 6 : Effectuer l'analyse QTL avec TR_RATE
  print("\n🔹 Analyse QTL avec TR_RATE...")

# Trouver l'index de la colonne TR_RATE
tr_rate_col <- which(colnames(cr_obj$pheno) == "BLUE_Slope")

# Pré-calculer les probabilités génotypiques
cr_obj <- calc.genoprob(cr_obj)

# Puis exécuter scanone
SIM <- scanone(
  cross = cr_obj,
  pheno.col = tr_rate_col,
  method = "hk"
)
# Visualiser les résultats
plot(SIM, col = "blue", main = "QTL Analysis for Slope")


#VOIR AVEC LES SEUILS


# Résumé de l'analyse QTL (inclut les seuils)
summary_SIM <- summary(SIM)

# Extraire les seuils
lod_standard <- summary_SIM$lod.thresholds["standard"]  # Seuil standard (ex: 5% LOD)
lod_bonferroni <- summary_SIM$lod.thresholds["bonferroni"]  # Seuil de Bonferroni
lod_permutation <- summary_SIM$lod.thresholds["permutation"]  # Seuil par permutation (si disponible)

print(paste("Seuil standard (5%) :", lod_standard))
print(paste("Seuil Bonferroni :", lod_bonferroni))


# Tracer le graphique avec les seuils
plot(
  SIM,
  col = "blue",
  main = "QTL Analysis for Slope",
  lodthresholds = c(
    "standard" = lod_standard,      # Seuil standard (5%)
    "bonferroni" = lod_bonferroni   # Seuil de Bonferroni
  )
)


# Tracer le graphique sans seuils
plot(SIM, col = "blue", main = "QTL Analysis for Slope")

# Ajouter les seuils manuellement
abline(
  h = lod_standard,
  col = "red",
  lty = 2,  # Ligne pointillée
  lwd = 2
)
abline(
  h = lod_bonferroni,
  col = "green",
  lty = 2,  # Ligne pointillée
  lwd = 2
)

# Ajouter une légende
legend(
  "topright",
  legend = c("Seuil standard (5%)", "Seuil Bonferroni"),
  col = c("red", "green"),
  lty = 2,
  lwd = 2,
  bty = "n"
)



#============================================================================================================================
#meme script que marcel en bas
#============================================================================================================================


set.seed(123)  # pour reproductibilité

perm <- scanone(cross = cr_obj,
                pheno.col = "BLUE_Slope",
                method = "hk",
                n.perm = 500,
                verbose = TRUE)



# voir les seuils
thresholds <- summary(perm, alpha = c(0.05, 0.01))
thresholds

#enrégistrer les 2 seuils calculés
write.table(thresholds, file = "thresholds_slope_BC12.txt", sep = "\t", quote = FALSE)
saveRDS(thresholds, file = "thresholds_slope_BC12.rds")

#Télécharger les seuils
thresholds <- readRDS("thresholds_slope_BC12.rds")

#LOD thresholds (1000 permutations)
#lod
#5% 4.36
#1% 5.04




# perform the QTL scan
SIM <- scanone(cross = cr_obj, pheno.col = "BLUE_Slope", method = "hk")
plot(SIM, col = "blue")
abline(h = thresholds[1, "lod"], col = "red", lty = 2, lwd = 2)


#si je veux mettre les 2 seuils sur le manhattan plot
plot(SIM, col = "blue")
abline(h = thresholds, col = "red", lty = 2, lwd = 2)
legend("topright", legend = paste("Seuil =", round(thresholds, 2)),
       col = "red", lty = 2, lwd = 2)

CIM <- cim(cross = cr_obj,
           pheno.col = "BLUE_Slope",
           method = "hk",
           n.marcovar = 5,
           window = 10)

plot(CIM, col = "blue")
abline(h = thresholds, col = "red", lty = 2, lwd = 2)
legend("topright", legend = paste("Seuil =", round(thresholds, 2)),
       col = "red", lty = 2, lwd = 2)


























#Calculer le seuil par permutation

set.seed(123)  # pour reproductibilité

perm <- scanone(cross = cr_obj_i,
                pheno.col = "RhizoSheathSize_C_ratio",
                method = "hk",
                n.perm = 1000,
                verbose = TRUE)

# voir les seuils
thresholds <- summary(perm, alpha = c(0.05, 0.01))
thresholds

#enrégistrer les 2 seuils calculés
write.table(thresholds, file = "thresholds_effet_fixe.txt", sep = "\t", quote = FALSE)
saveRDS(thresholds, file = "thresholds_effet_fixe.rds")

#Télécharger les seuils
thresholds <- readRDS("thresholds_effet_fixe.rds")

#LOD thresholds (1000 permutations)
#lod
#5% 4.13
#1% 4.95


#Identification des qtls significatifs
#QTLs <- summary(CIM, threshold = 3.80)
QTLs_CIM <- summary(CIM, threshold = thresholds[1, "lod"])
#thresholds["5%", "lod"]

# Vérifier le QTL sélectionné
#print(QTLs)
QTLs_CIM

find.marker(cr_obj_i, chr = 4, pos = 136)


QTL_final <- QTLs_CIM
colnames(QTL_final) <- c("Chr", "Position", "LOD")
QTL_final

#L'intervalle de confiance
lodint(CIM, chr = S4_66061444, drop = 1.5)


#version automatique pour tous les marqueurs
peaks <- QTLs_CIM

IC_list <- lapply(1:nrow(peaks), function(i) {
  lodint(CIM,
         chr = peaks$chr[i],
         drop = 1.5)
})

IC_list

#Pour avoir le R2
qtl_model <- makeqtl(cr_obj_i,
                     chr = c(6),
                     pos = c(69.23302),
                     what = "prob")


fit <- fitqtl(cr_obj_i,
              qtl = qtl_model,
              pheno.col = "RhizoSheathSize_C_ratio",
              method = "hk")

summary(fit)




#Colocalisation
IC_table <- data.frame(
  QTL = rownames(QTLs_CIM),
  chr = QTLs_CIM$chr,
  pos = QTLs_CIM$pos,
  stringsAsFactors = FALSE
)

IC_table$IC <- lapply(1:nrow(IC_table), function(i) {
  lodint(CIM,
         chr = IC_table$chr[i],
         drop = 1.5)
})


IC_bounds <- do.call(rbind, lapply(IC_table$IC, function(x) {
  data.frame(
    chr = unique(x$chr),
    start = min(x$pos),
    end = max(x$pos)
  )
}))

IC_bounds$QTL <- IC_table$QTL
IC_bounds


overlap <- function(a_start, a_end, b_start, b_end) {
  return(!(a_end < b_start | b_end < a_start))
}

n <- nrow(IC_bounds)
coloc_mat <- matrix(FALSE, n, n)

for (i in 1:n) {
  for (j in 1:n) {
    if (IC_bounds$chr[i] == IC_bounds$chr[j]) {
      coloc_mat[i, j] <- overlap(IC_bounds$start[i], IC_bounds$end[i],
                                 IC_bounds$start[j], IC_bounds$end[j])
    }
  }
}

coloc_mat


#Colocalisation avec gapit
IC_bounds


effectplot(cr_obj_i,
           pheno.col = "RhizoSheathSize_C_ratio",
           mname1 = find.marker(cr_obj_i, chr = 6, pos = 69.23302))



#Message d'avis :
#Dans effectplot(cr_obj_i, pheno.col = 2, mname1 = find.marker(cr_obj_i,  :
#   -Running sim.geno.

cr_obj_i <- sim.geno(cr_obj_i, n.draws = 100)

effectplot(cr_obj_i,
           pheno.col = "RhizoSheathSize_C_ratio",
           mname1 = find.marker(cr_obj_i, chr = 6, pos = 69.23302))

#effet QTL
effectplot(cr_obj_i, pheno.col="RhizoSheathSize_C_ratio", mname1 = "S6_49256435" )      
plotPXG(cr_obj_i,marker ="S6_49256435",pheno.col="RhizoSheathSize_C_ratio")

#je veux déterminer l'effet qu QTL comme dans mppR
marker <- find.marker(cr_obj_i, chr = 6, pos = 69.23302)

g <- pull.geno(cr_obj_i)[, marker]

means <- tapply(cr_obj_i$pheno$RhizoSheathSize_C_ratio, g, mean)
sdv   <- tapply(cr_obj_i$pheno$RhizoSheathSize_C_ratio, g, sd)
n     <- table(g)

se <- sdv / sqrt(n)

res <- data.frame(
  Genotype = names(means),
  Effect = means - means["1"],   # référence = parent 1 (Grinkan)
  StdErr = se
)

res

#ici, n est convertit en numérique
n <- as.numeric(table(g))

se <- sdv / sqrt(n)

res <- data.frame(
  Genotype = names(means),
  Mean = means,
  Effect = means - means["1"],
  SE = se,
  N = n
)

res

g <- pull.geno(cr_obj_i)[, marker]

model <- lm(RhizoSheathSize_C_ratio ~ as.factor(g), data = cr_obj_i$pheno)

summary(model)
































































# perform the QTL scan
SIM <- scanone(cross = cr_obj_i, pheno.col = "RhizoSheathSize_C_ratio", method = "hk")
plot(SIM, col = "blue")
abline(h = thresholds[1, "lod"], col = "red", lty = 2, lwd = 2)


#si je veux mettre les 2 seuils sur le manhattan plot
plot(SIM, col = "blue")
abline(h = thresholds, col = "red", lty = 2, lwd = 2)
legend("topright", legend = paste("Seuil =", round(thresholds, 2)),
       col = "red", lty = 2, lwd = 2)


#Identification des qtls significatifs
QTLs <- summary(SIM, threshold = 4.58)

# Vérifier le QTL sélectionné
print(QTLs)


CIM <- cim(cross = cr_obj_i,
           pheno.col = "RhizoSheathSize_C_ratio",
           method = "hk",
           n.marcovar = 5,
           window = 10)

plot(CIM, col = "blue")
abline(h = thresholds, col = "red", lty = 2, lwd = 2)
legend("topright", legend = paste("Seuil =", round(thresholds, 2)),
       col = "red", lty = 2, lwd = 2)



#Identification des qtls significatifs
#QTLs <- summary(CIM, threshold = 3.80)
QTLs_CIM <- summary(CIM, threshold = thresholds[1, "lod"])
#thresholds["5%", "lod"]

# Vérifier le QTL sélectionné
#print(QTLs)
QTLs_CIM

find.marker(cr_obj_i, chr = 4, pos = 136)


QTL_final <- QTLs_CIM
colnames(QTL_final) <- c("Chr", "Position", "LOD")
QTL_final

#L'intervalle de confiance
lodint(CIM, chr = 6, drop = 1.5)


#version automatique pour tous les marqueurs
peaks <- QTLs_CIM

IC_list <- lapply(1:nrow(peaks), function(i) {
  lodint(CIM,
         chr = peaks$chr[i],
         drop = 1.5)
})

IC_list

#Pour avoir le R2
qtl_model <- makeqtl(cr_obj_i,
                     chr = c(6),
                     pos = c(69.23302),
                     what = "prob")


fit <- fitqtl(cr_obj_i,
              qtl = qtl_model,
              pheno.col = "RhizoSheathSize_C_ratio",
              method = "hk")

summary(fit)




#Colocalisation
IC_table <- data.frame(
  QTL = rownames(QTLs_CIM),
  chr = QTLs_CIM$chr,
  pos = QTLs_CIM$pos,
  stringsAsFactors = FALSE
)

IC_table$IC <- lapply(1:nrow(IC_table), function(i) {
  lodint(CIM,
         chr = IC_table$chr[i],
         drop = 1.5)
})


IC_bounds <- do.call(rbind, lapply(IC_table$IC, function(x) {
  data.frame(
    chr = unique(x$chr),
    start = min(x$pos),
    end = max(x$pos)
  )
}))

IC_bounds$QTL <- IC_table$QTL
IC_bounds


overlap <- function(a_start, a_end, b_start, b_end) {
  return(!(a_end < b_start | b_end < a_start))
}

n <- nrow(IC_bounds)
coloc_mat <- matrix(FALSE, n, n)

for (i in 1:n) {
  for (j in 1:n) {
    if (IC_bounds$chr[i] == IC_bounds$chr[j]) {
      coloc_mat[i, j] <- overlap(IC_bounds$start[i], IC_bounds$end[i],
                                 IC_bounds$start[j], IC_bounds$end[j])
    }
  }
}

coloc_mat


#Colocalisation avec gapit
IC_bounds


effectplot(cr_obj_i,
           pheno.col = "RhizoSheathSize_C_ratio",
           mname1 = find.marker(cr_obj_i, chr = 6, pos = 69.23302))



#Message d'avis :
#Dans effectplot(cr_obj_i, pheno.col = 2, mname1 = find.marker(cr_obj_i,  :
#   -Running sim.geno.

cr_obj_i <- sim.geno(cr_obj_i, n.draws = 100)

effectplot(cr_obj_i,
           pheno.col = "RhizoSheathSize_C_ratio",
           mname1 = find.marker(cr_obj_i, chr = 6, pos = 69.23302))

#effet QTL
effectplot(cr_obj_i, pheno.col="RhizoSheathSize_C_ratio", mname1 = "S6_49256435" )      
plotPXG(cr_obj_i,marker ="S6_49256435",pheno.col="RhizoSheathSize_C_ratio")

#je veux déterminer l'effet qu QTL comme dans mppR
marker <- find.marker(cr_obj_i, chr = 6, pos = 69.23302)

g <- pull.geno(cr_obj_i)[, marker]

means <- tapply(cr_obj_i$pheno$RhizoSheathSize_C_ratio, g, mean)
sdv   <- tapply(cr_obj_i$pheno$RhizoSheathSize_C_ratio, g, sd)
n     <- table(g)

se <- sdv / sqrt(n)

res <- data.frame(
  Genotype = names(means),
  Effect = means - means["1"],   # référence = parent 1 (Grinkan)
  StdErr = se
)

res

#ici, n est convertit en numérique
n <- as.numeric(table(g))

se <- sdv / sqrt(n)

res <- data.frame(
  Genotype = names(means),
  Mean = means,
  Effect = means - means["1"],
  SE = se,
  N = n
)

res

g <- pull.geno(cr_obj_i)[, marker]

model <- lm(RhizoSheathSize_C_ratio ~ as.factor(g), data = cr_obj_i$pheno)

summary(model)



















