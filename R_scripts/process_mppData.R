##########################
# Process mppData object #
##########################

# library ----

library(mppR)
library(qtl)
# library(ASMap)
library(dplyr)
library(data.table)

# ad-hoc function ----
source("./functions/IBD_mppData_BCNAM.R")

# compile the ABH matrix ----

# Get parents
parents <- read.table(file = 'D:/BCNAM_data/05_Global_matrices/All_parents_IUPAC.txt', header = TRUE)

n_mk <- nrow(parents)

het_sc_lk <- c(0, 1, 1, 1, NA, 2, 1, 2, 2)
names(het_sc_lk) <- c('00', '01', '02', '10', '11', '12', '20', '21', '22')

# Global par_per_cross
# pop_id <- c('GR', 'KK', 'Lata')
# p <- 1
# load(file = paste0('data/par_per_cross/par_per_cross_', pop_id[p], '.RData'))

# Iterate over the crosses self-defined
n_cr <- 2

folder_raw_data <- "D:/BCNAM_data/04_2_Final_RQTL_input/"
pop_id <- 'BC12-31'

par_per_cross <- data.frame(cross = c("BC12", "BC31"),
                            Par1 = rep("GR", n_cr),
                            Par2 = c("Lata3", "B35"),
                            p1_id = rep("v12", n_cr),
                            p2_id = c("X1097", "v33"))

# create some space
geno_ABH <- matrix(NA, nrow = n_cr*150, ncol = n_mk)
geno_012 <- matrix(NA, nrow = n_cr*150, ncol = n_mk)
allele_freq <- matrix(NA, n_cr, 3)
colnames(allele_freq) <- c('A', 'B', 'H')
rownames(allele_freq) <- par_per_cross$cross
geno_id_vec <- c()
  
f_row <- 1
  
# setwd(folder_raw_data)
  
# Iterate over the crosses of the population
t1 <- Sys.time()
  
for(c in 1:n_cr){
    
    cr <- par_per_cross$cross[c]
    RP <- par_per_cross$p1_id[c]
    DP <- par_per_cross$p2_id[c]
    
    # Get ABH matrix
    # RQTL_fi <- read.table(file = paste0(cr, "_final_RQTL_input.txt"), sep = ",")
    RQTL_fi <- fread(file = paste0(folder_raw_data, cr, "_final_RQTL_input.txt"),
                     sep = ",", header = FALSE)
    RQTL_fi <- data.frame(RQTL_fi)
    Y <- RQTL_fi[4:nrow(RQTL_fi), 4:ncol(RQTL_fi)]
    Y <- t(Y)
    
    # Fill the remaining missing values
    for(j in 1:ncol(Y)){
      
      Y_miss <- Y[, j]
      if(any(is.na(Y_miss))){
        
        miss_pos <- which(is.na(Y_miss))
        Y_miss[miss_pos] <- Y_miss[miss_pos-1]
        Y[, j] <- Y_miss
        
      }
      
    }
    
    # Check if there are missing values
    if(sum(is.na(c(as.matrix(Y)))) !=0){stop('missing values')}
    
    # fill the allele frequency information
    allele_freq[c, ] <- round(100*table(c(Y))/(nrow(Y) * ncol(Y)), 1)
    
    Site <- unlist(RQTL_fi[1, 4:ncol(RQTL_fi)])
    geno_id <- RQTL_fi[4:nrow(RQTL_fi), 1]
    Y <- data.frame(Site, Y)
    colnames(Y) <- c('Site', geno_id)
    
    BCXX_complete <- merge(x=parents[, c("Site", "alleles", "chrom", "pos", RP, DP)],
                           y = Y, by = "Site", all.x = TRUE)
    
    # Ordonne la matrice
    BCXX_complete <- BCXX_complete[order(BCXX_complete$chrom, BCXX_complete$pos), ]
    
    # Parents conversion
    geno_par <- data.frame(A=BCXX_complete[, 5], B=BCXX_complete[, 6])
    ref_allele <- substr(BCXX_complete$alleles, 1, 1)
    alt_allele <- substr(BCXX_complete$alleles, 3, 3)
    
    # Convert RP (parent A) into ref or alt
    RP_allele <- rep(NA, n_mk)
    
    test_ref <- geno_par$A == ref_allele
    test_alt <- geno_par$A == alt_allele
    test_het <- !test_ref & !test_alt
    
    RP_allele[test_ref] <- 0
    RP_allele[test_alt] <- 2
    RP_allele[test_het] <- 1
    
    # Convert DP (parent B) into ref or alt
    DP_allele <- rep(NA, n_mk)
    
    test_ref <- geno_par$B == ref_allele
    test_alt <- geno_par$B == alt_allele
    test_het <- !test_ref & !test_alt
    
    DP_allele[test_ref] <- 0
    DP_allele[test_alt] <- 2
    DP_allele[test_het] <- 1
    
    het_allele <- as.character(paste0(RP_allele, DP_allele))
    
    par_allele <- cbind(RP_allele, DP_allele, het_sc_lk[het_allele])
    
    colnames(par_allele) <- c("A", "B", "H")
    
    BCXX_complete[, 5:6] <- par_allele[, 1:2]
    
    matrix_012 <- matrix(NA, n_mk, ncol(BCXX_complete) - 6)
    
    # transform the genotype matrix into 012
    for(i in 1:n_mk){
      
      mk_lk <- par_allele[i, ]
      mk_sc <- as.character(BCXX_complete[i, 7:ncol(BCXX_complete)])
      matrix_012[i, ] <- mk_lk[mk_sc]
      
    }
    
    # fill the reference matrix
    end_ind <- f_row + (ncol(Y) - 1)
    
    geno_ABH[f_row:(end_ind - 1), ] <- t(Y[, -1])
    geno_012[f_row:(end_ind - 1), ] <- t(matrix_012)
    geno_id_vec <- c(geno_id_vec, geno_id) 
    
    f_row <- end_ind
    cat("End: ", par_per_cross$cross[c])
    
  }
  
t2 <- Sys.time()
print(t2-t1)
  
# Save the two global matrices as well as the allele frequencies
geno_ABH <- geno_ABH[!is.na(geno_ABH[, 1]), ]
geno_012 <- geno_012[!is.na(geno_012[, 1]), ]
  
colnames(geno_ABH) <- colnames(geno_012) <- parents$Site
rownames(geno_ABH) <- rownames(geno_012) <- geno_id_vec
  
setwd('C:/Users/vince/OneDrive/Documents/WD/ICRISAT/BCNAM/data/genotype')
  
save(geno_ABH, file = file.path('data/geno', 'geno_ABH',
                                paste0('geno_ABH_', pop_id,'.RData')))

save(geno_012, file = file.path('data/geno', 'geno_012',
                                paste0('geno_012_', pop_id,'.RData')))

# save(allele_freq, file = paste0('allele_freq_', pop_id[p], '.RData'))
  


# create the mppData object ----

# par_per_cross
pop_id <- 'BC12-31'
n_cr <- 2
par_per_cross <- data.frame(cross = c("BC12", "BC31"),
                            Par1 = rep("GR", n_cr),
                            Par2 = c("Lata3", "B35"),
                            p1_id = rep("v12", n_cr),
                            p2_id = c("X1097", "v33"))

par_per_cross <- as.matrix(par_per_cross)

# map
load(file = 'data/map/Global_map.Rdata')
map <- map[, -4]
n_mk <- nrow(map)

# geno offsping (ABH)
load(file = file.path('data/geno', 'geno_ABH',
                      paste0('geno_ABH_', pop_id,'.RData')))

rownames(geno_ABH) <- toupper(rownames(geno_ABH))

# check that the list of markers is the same
stopifnot(identical(map$mk.names, colnames(geno_ABH)))

# geno parents
RP <- rep('A', n_mk)
DP <- matrix('B', nrow = nrow(par_per_cross), ncol = n_mk)
geno_par <- rbind(RP, DP)
colnames(geno_par) <- map$mk.names
rownames(geno_par)<- c(par_per_cross[1, 2], par_per_cross[, 3])

# pheno data: load here phenotype data
pheno <- as.matrix(data.frame(t1 = rnorm(nrow(geno_ABH))))
rownames(pheno) <- rownames(geno_ABH)

# cross indicator
cross_ind <- substr(rownames(pheno), 1, 4)

stopifnot(identical(rownames(geno_ABH), rownames(pheno)))

# create the mppData object
mppData <- create.mppData(geno.off = geno_ABH, geno.par = geno_par, map = map,
                          pheno = pheno, cross.ind = cross_ind,
                          par.per.cross = par_per_cross)

# skip QC of the genotype data
mppData$status <- 'QC'
rm(geno_ABH)

# add geno.id and ped.mat
geno.id <- rownames(mppData$pheno)
mppData$geno.id <- geno.id

p1 <- mppData$par.per.cross[, 2]
p2 <- mppData$par.per.cross[, 3]

cross.ind <- mppData$cross.ind

names(p2) <- names(p1) <- mppData$par.per.cross[, 1]

ped.mat <- data.frame(rep("offspring", length(geno.id)), geno.id,
                      p1[cross.ind], p2[cross.ind], stringsAsFactors = FALSE)
colnames(ped.mat) <- c("type" ,"genotypes", "parent1", "parent2")

mppData$ped.mat <- ped.mat

# add manually the 012 data

# geno offsping (012)
load(file = file.path('data/geno', 'geno_012',
                      paste0('geno_012_', pop_id,'.RData')))
rownames(geno_012) <- toupper(rownames(geno_012))

geno_012 <- geno_012[rownames(mppData$geno.off), ]

stopifnot(identical(rownames(geno_012), rownames(mppData$pheno)))

mppData$geno.IBS <- geno_012
mppData$status <- 'IBS'
rm(geno_012)

# IBD mppData
mppData <- IBD.mppData_BCNAM(mppData = mppData)

rm(geno_par)
rm(pheno)

# count the number of cross-over
# CO_count <- CO_stat(mppData)
# setwd('C:/Users/vince/OneDrive/Documents/WD/ICRISAT/BCNAM/Results/Genetic_analysis/CO_count')
# save(CO_count, file = paste0('CO_count_', pop_id[i], '.RData'))

# save the mppData object global

save(mppData, file = file.path('data/mppData',
                               paste0('mppData_', pop_id,'.RData')))