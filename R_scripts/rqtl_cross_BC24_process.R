##################################
# mppData modification pour rqtl #
##################################

# library
library(mppR)
library(qtl)
library(readxl)

# data

# pheno <- read_xlsx(path = "data/pheno/Data_rhizosheath_BC24.xlsx")
# load(file = "data/map/map_BC24_filtrage.RData")
# # load(file = "data/mppData/mppData.RData")
# # 
# # # check the genotype
# # geno.id <- mppData$geno.id
# # 
# # all(pheno$Lines %in% geno.id)
# # sum(!(pheno$Lines %in% geno.id))
# # 
# # mppData$par.per.cross
# 
# # the mppData object does not contain BC24 cross
# 
# # need to start from the initial mppData object
# load(file = "data/mppData/raw_data/GR/mppData.RData")
# 
# # check is the cross BC24 is present
# mppData$par.per.cross
# 
# # ok
# 
# # check the genotype
# geno.id <- mppData$geno.id
# 
# all(pheno$Lines %in% geno.id)
# sum(!(pheno$Lines %in% geno.id))
# 
# # not all genotype have genotypic information
# 
# # determine common genotype
# geno_com <- intersect(pheno$Lines, mppData$geno.id)
# 
# # filter the mppData object to keep only the genotype from BC24
# mppData <- subset(mppData, gen.list = geno_com)
# 
# # filtrer le mppData objet avec les marqueurs présents pour la GWAS
# mppData <- subset(mppData, mk.list = map2$SNP)
# 
# # récupérer l'objet cross pour faire l'analyse avec r/qtl
# geno_rqtl <- mppData$geno.IBD
# 
# # modification de l'objet pheno
# pheno_org <- geno_rqtl$pheno
# 
# pheno_red <- pheno[match(geno_com, pheno$Lines), ]
# pheno_rqtl <- pheno_org[1:nrow(pheno_red), 1, drop = FALSE]
# pheno_rqtl$SB1_FLAG <- pheno_red$RAS.RT
# 
# colnames(pheno_rqtl)[1] <- colnames(pheno_red)[2]
# geno_rqtl$pheno <- pheno_rqtl
# 
# # analyses : SIM, CIM, ...
# sim <- scanone(cross = geno_rqtl)
# 
# # Does not work. Alternative strategy.
# 
# # 1) subset the marker data from the mppData
# pheno <- read_xlsx(path = "data/pheno/Data_rhizosheath_BC24.xlsx")
# map <- load(file = "data/map/map_BC24_filtrage.RData")
# load(file = "data/mppData/raw_data/GR/mppData.RData")
# 
# mppData <- subset(mppData, mk.list = map2$SNP)
# 
# # 2) get the r/qtl object
# cross_obj <- mppData$geno.IBD
# 
# 
# # 3) subset with the common genotypes
# 
# # determine common genotype
# geno_com <- intersect(pheno$Lines, mppData$geno.id)
# geno_sel <- mppData$geno.id %in% geno_com 
# 
# # subset the cross object
# cross_obj <- cross_obj[, geno_sel]
# 
# # modify the phenotype values
# pheno_red <- pheno[match(geno_com, pheno$Lines), ]
# cross_obj$pheno[, 1] <- pheno_red[, 2]
# 
# # analyses : SIM, CIM, ...
# sim <- scanone(cross = cross_obj)

# try from the script of the individual analyses

pheno <- read_xlsx(path = "data/pheno/Data_rhizosheath_BC24.xlsx")
pheno <- as.data.frame(pheno)
load(file = "data/map/map_BC24_filtrage.RData")

rownames(pheno) <- pheno$Lines

geno_pheno <- pheno$Lines
geno_pheno <- geno_pheno[!is.na(geno_pheno)]

# mppData
load(file = "data/mppData/raw_data/GR/mppData.RData")

# extract the rqtl object from the mppData object
cr_obj <- mppData$geno.IBD

# subset of the reduced list of markers used in GWAS
mk_to_remove <- mppData$map$mk.names[!(mppData$map$mk.names %in% map2$SNP)]
cr_obj <- drop.markers(cross = cr_obj, markers = mk_to_remove)
cr_obj <- calc.genoprob(cross = cr_obj)

# possibilité de sauver le cross object ici pour ne pas devoir le
# reproduire à chaque fois.

# list of genotypes in the large mppData object (all Grinkan)
geno_id_ref <- mppData$geno.id

# common genotype between phenotype and genotypes
geno_common <- base::intersect(mppData$geno.id, geno_pheno)

# subset the cross object using the common genotypes
cr_obj_i <- subset(x = cr_obj, ind = geno_id_ref %in% geno_common)

# subset the phenotype data using the common genotypes
pheno_i <- pheno[geno_common, ]
rownames(pheno_i) <- NULL

# add the phenotype to the cross object
cr_obj_i$pheno <- pheno_i

# perform the QTL scan
SIM <- scanone(cross = cr_obj_i, pheno.col = 2, method = "hk")
plot(SIM, col = "blue")

QTLs <- summary(SIM, threshold = 2.5)
