###############################
# mppData object construction #
###############################

library(mppR)

# determine from which recurrent parent the cross come from ----

# list of genotpes
# load(file = "output/geno/liste_Slope_genotypes.RData")
# 
# geno_id_Slope <- genotypes_BLUE
# 
# load(file = "output/geno/liste_TR_RATE_24H_genotypes.RData")
# 
# geno_id_TR_RATE_24H <- genotypes_BLUE
# 
# load(file = "output/geno/liste_genotypes_en_pots_2026.RData")

geno_id <- unique(union(union(as.character(geno_id_Slope),
                              as.character(geno_id_TR_RATE_24H)), as.character(all_genotypes)))

# récupérer liste genotype BC12 et BC31



cross_ind <- substr(x = geno_id, 1, 4)

# par_per_cross
load(file = "data/par_per_cross/par_per_cross_KK.RData")
ppc_KK <- par_per_cross

load(file = "data/par_per_cross/par_per_cross_GR.RData")
ppc_GR <- par_per_cross


#GR x Lata_3 etant absent, je l'ajoute 
# Ajouter la nouvelle croissement GR x Lata_3 avec ids parents
nouvelle_population <- data.frame(
  cross = "BC12",
  Par1  = "GR",
  Par2  = "Lata_3",
  p1_id = "v12",
  p2_id = "v14",
  stringsAsFactors = FALSE
)
ppc_GR <- rbind(ppc_GR, nouvelle_population)

par_per_cross <- rbind(ppc_KK, ppc_GR)

# Restreindre aux crosses présentes dans tes génotypes
ppc_data <- par_per_cross[par_per_cross$cross %in% cross_ind, ]

unique(ppc_data$Par1)

# Only data from the Grinkan sub-population.

# Subset from the Grinkan global mppData ----

# récupérer séparément les données génétiques de BC12 et BC31
# et les mettre ensemble.

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
# save(mppData, file = "data/mppData/mppData.RData")
# Vincent m'a demandé de ne pas exécuter la ligne 75

