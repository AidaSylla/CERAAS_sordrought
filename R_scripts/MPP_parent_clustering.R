#########################################
# clustering of parent using clusthaplo #
#########################################

library(mppR)
library(clusthaplo)

load(file = "data/mppData/mppData.RData")

mppData$par.per.cross

haplo.map <- mppData$map[, -3]
consensus.map <- mppData$map[, -3]
map <- mppData$map
parents <- mppData$parents

code_parents <- c("v3", "v12", "v18", "v25")

# Introduce the genotype of the parents

# example
# data(USNAM_geno)
# geno.par <- USNAM_geno[1:6, ]

# parent data from BCNAM
load(file = "D:/Mes Donnees/WD/BCNAM/data/genotype/all_parents_geno.RData")
geno_par <- data[, code_parents]
colnames(geno_par) <- c("CSM388", "Grinkan", "IS23540", "White_Kaur")

# marker.data <- t(mppData$geno.par.clu)
marker.data <- geno_par
rownames(marker.data) <- map$mk.names

# clusthaplo execution

chr.fact <- factor(x = map[, 2], levels = unique(map[, 2]))

step.size <- max(tapply(X = map[, 4], INDEX = chr.fact, FUN = max)) + 100 

# function default arguments
w1 = "kernel.exp"
w2 = "kernel.unif"
window = 25
K = 10
simulation.type = "equi"
simulation.Ng = 50
simulation.Nrep = 2
threshold.quantile = 95
plot = TRUE
plot.loc = getwd()

p_clu <- parent_cluster(haplo.map = haplo.map,
                        consensus.map = consensus.map,
                        marker.data = marker.data,
                        na.strings = NA,
                        w1 = w1, w2 = w2,
                        step.size = step.size,
                        window = 10, K = K,
                        simulation.type = simulation.type,
                        simulation.Ng = simulation.Ng,
                        simulation.Nrep = simulation.Nrep,
                        threshold.quantile = threshold.quantile,
                        plot = plot,
                        plot.loc = plot.loc)

p_clu_w10 <- p_clu$par.clu

p_clu <- parent_cluster(haplo.map = haplo.map,
                        consensus.map = consensus.map,
                        marker.data = marker.data,
                        na.strings = NA,
                        w1 = w1, w2 = w2,
                        step.size = step.size,
                        window = 25, K = K,
                        simulation.type = simulation.type,
                        simulation.Ng = simulation.Ng,
                        simulation.Nrep = 2,
                        threshold.quantile = threshold.quantile,
                        plot = plot,
                        plot.loc = plot.loc)

p_clu_w25 <- p_clu$par.clu

# nombre d'allèles ancestraux moyen doit être compris entre 1 et n_parents
# [1, 4]
# windows = 10 => 1.03
# windows = 25 => 1.00225
# windows = 5 => ?
# windows = 2 => ?

p_clu$av.cl

# nouveaux modèles après avoir fait le clustering the parents

mppData$par.clu <- p_clu_w10
MPP_SIM_anc <- mpp_SIM(mppData = mppData, Q.eff = "anc")
MPP_mQTL <- MQE_forward(mppData = mppData, Q.eff = c("par", "anc"),
                        threshold = 3)