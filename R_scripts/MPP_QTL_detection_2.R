######################
# Script MPP analyse #
######################

# library ----
library(mppR)

# load data ----
load(file = "data/mppData/mppData.RData")
load(file = "output/parent_cluster/par_clu.RData")
mppData$par.clu <- par_clu
mppData$status <- "complete"

# subset des données ----
# Réduire la taille des données pour tester plus rapidement

# selection 1000 le long du génome
#mk_sel <- sort(sample(x = 1:nrow(mppData$map), size = 1000))
#mppData <- subset(mppData, mk.list = mk_sel)

# trait
traits <- colnames(mppData$pheno)
n_traits <- length(traits)

# modèles: par, anc, MQE (par, anc)
models <- c("par", "anc", "MQE")
n_models <- length(models)

for (i in 1:n_traits){
  
  for(j in 1:n_models){
    
    if(models[j] == "MQE"){
      m <- MQE_proc(pop.name = "BCNAM", trait.name = traits[i],
                    mppData = mppData, trait = i, Q.eff = c("par", "anc"),
                    output.loc = "results/MPP_QTL_analysis")
    } else{
      m <- mpp_proc(pop.name = "BCNAM", trait.name = traits[i],
                    mppData = mppData, trait = i, Q.eff = models[j],
                    output.loc = "results/MPP_QTL_analysis")
    }
    
  }
  
}



















#EN BAS C'EST MOI


































######################
# Script MPP analyse #
######################

# library ----
library(mppR)

# load data ----
load(file = "data/mppData/mppData.RData")
load(file = "output_Aida_all_snps/parent_cluster/par_clu.RData")
mppData$par.clu <- par_clu
mppData$status <- "complete"

# subset des données ----
# Réduire la taille des données pour tester plus rapidement

# selection 1000 le long du génome
#mk_sel <- sort(sample(x = 1:nrow(mppData$map), size = 1000))
#mppData <- subset(mppData, mk.list = mk_sel)

# trait
traits <- colnames(mppData$pheno)
n_traits <- length(traits)

# modèles: par, anc, MQE (par, anc)
models <- c("par", "anc", "MQE")
n_models <- length(models)

for (i in 1:n_traits){
  
  for(j in 1:n_models){
    
    if(models[j] == "MQE"){
      m <- MQE_proc(pop.name = "BCNAM", trait.name = traits[i],
                    mppData = mppData, trait = i, Q.eff = c("par", "anc"),
                    output.loc = "output_Aida/results/MPP_QTL_analysis")
    } else{
      m <- mpp_proc(pop.name = "BCNAM", trait.name = traits[i],
                    mppData = mppData, trait = i, Q.eff = models[j],
                    output.loc = "output_Aida/results/MPP_QTL_analysis")
    }
    
  }
  
}

