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
