#########################################
# clustering of parent using clusthaplo #
#########################################
#pour installer clusthaplo je vais dans packages puis installer puis package archive puis dans le projet puis function puis library

library(mppR)
library(clusthaplo)

#load(file = "data/mppData/mppData.RData")
#load(file = "data/mppData/raw_data/GR/mppData.RData")
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
#load(file = "D:/Mes Donnees/WD/BCNAM/data/genotype/all_parents_geno.RData")
load(file = "C:/Users/2025an002/Desktop/Sordrought_2025/CERAAS_sordrought/data/geno/all_parents_geno.RData")

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
p_clu$av.cl

p_clu <- parent_cluster(haplo.map = haplo.map,
                        consensus.map = consensus.map,
                        marker.data = marker.data,
                        na.strings = NA,
                        w1 = w1, w2 = w2,
                        step.size = step.size,
                        window = 5, K = K,
                        simulation.type = simulation.type,
                        simulation.Ng = simulation.Ng,
                        simulation.Nrep = 2,
                        threshold.quantile = threshold.quantile,
                        plot = plot,
                        plot.loc = plot.loc)

p_clu_w5 <- p_clu$par.clu
p_clu$av.cl

p_clu <- parent_cluster(haplo.map = haplo.map,
                        consensus.map = consensus.map,
                        marker.data = marker.data,
                        na.strings = NA,
                        w1 = w1, w2 = w2,
                        step.size = step.size,
                        window = 2, K = K,
                        simulation.type = simulation.type,
                        simulation.Ng = simulation.Ng,
                        simulation.Nrep = 2,
                        threshold.quantile = threshold.quantile,
                        plot = plot,
                        plot.loc = plot.loc)

p_clu_w2 <- p_clu$par.clu
p_clu$av.cl

p_clu <- parent_cluster(haplo.map = haplo.map,
                        consensus.map = consensus.map,
                        marker.data = marker.data,
                        na.strings = NA,
                        w1 = w1, w2 = w2,
                        step.size = step.size,
                        window = 1, K = K,
                        simulation.type = simulation.type,
                        simulation.Ng = simulation.Ng,
                        simulation.Nrep = 2,
                        threshold.quantile = threshold.quantile,
                        plot = plot,
                        plot.loc = plot.loc)

p_clu_w1 <- p_clu$par.clu
p_clu$av.cl

#source("C:/Users/2025an002/Desktop/Sordrought_2025/CERAAS_sordrought/functions/fct_parent_cluster.R", echo = TRUE)


# nombre d'allèles ancestraux moyen doit être compris entre 1 et n_parents
# [1, 4]
# windows = 25 => 1.00225
# windows = 10 => 1.045048
# windows = 5 => 1.138985
# windows = 2 => 1.474168
# windows = 1 => 1.830032

# nombre d'allèles ancestraux moyen doit être compris entre 1 et n_parents que moi j'ai trouvé
# [1, 4]
# windows = 25 => 
# windows = 10 => 1.042526
# windows = 5 => 1.139587
# windows = 2 => 1.468387
# windows = 1 => 1.826967
p_clu$av.cl

# nouveaux modèles après avoir fait le clustering the parents

mppData$par.clu <- p_clu_w2
MPP_SIM_anc <- mpp_SIM(mppData = mppData, Q.eff = "anc")
MPP_mQTL <- MQE_forward(mppData = mppData, Q.eff = c("par", "anc"),
                        threshold = 3)

# Selection du clustering parental: On choisit le clustering parental donnée par
# win = 1. Sur la base d'un argument plutot statistique afin d'avoir suffisamment
# de groupement d'allèles parentaux le long du génome. Avec win = 1, on a en moyenne
# 1.83 allèles parentaux sur 4. Groupement moyen.

par_clu <- p_clu_w1

save(par_clu, file = "output/parent_cluster/par_clu.RData")
save(par_clu, file = "output_Aida/parent_cluster/par_clu.RData")




p_clu_w10 <- p_clu$par.clu
p_clu$av.cl
[1] 1.041478
p_clu_w5 <- p_clu$par.clu
> p_clu$av.cl
[1] 1.139781
p_clu_w2 <- p_clu$par.clu
> p_clu$av.cl
[1] 1.483519
p_clu_w1 <- p_clu$par.clu
> p_clu$av.cl
[1] 1.816743


