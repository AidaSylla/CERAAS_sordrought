######################
# Script MPP analyse #
######################

# library ----
library(mppR)

# load data ----
load(file = "data/mppData/mppData.RData")
load(file = "output_A/parent_cluster/par_clu.RData")

mppData$par.clu <- par_clu
mppData$status <- "complete"

# traits ----
traits <- colnames(mppData$pheno)
n_traits <- length(traits)

# modèles ----
models <- c("par", "anc", "MQE")
n_models <- length(models)

###############################
# Calcul des seuils permutation
###############################

perm_results <- list()

for(i in 1:n_traits){
  
  cat("Permutation pour :", traits[i], "\n")
  
  perm_results[[traits[i]]] <- mpp_perm(
    mppData = mppData,
    trait = i,
    N = 1000,      # nombre de permutations
    Q.eff = "par",      # ou "anc"
  )
  
}

###LE CACUL DU SEUIL AVEC LE SCRIPT EN BAS NE MARCHE PAS
for(i in 1:n_traits){
  
  cat("Permutation pour :", traits[i], "\n")
  
  perm_results[[traits[i]]] <- perm_test(
    mppData = mppData,
    trait = i,
    n.perm = 1000,      # nombre de permutations
    alpha = 0.05,       # seuil de significativité
    Q.eff = "par",      # ou "anc"
    plot = TRUE
  )
  
}
###################################
# Analyse QTL avec seuil permutation
###################################

for(i in 1:n_traits){
  
  # extraction du seuil
  thr <- perm_results[[traits[i]]]$threshold
  
  cat("Trait :", traits[i], 
      "- seuil permutation =", thr, "\n")
  
  for(j in 1:n_models){
    
    cat("Modèle :", models[j], "\n")
    
    if(models[j] == "MQE"){
      
      m <- MQE_proc(
        pop.name = "BCNAM",
        trait.name = traits[i],
        mppData = mppData,
        trait = i,
        Q.eff = c("par", "anc"),
        threshold = thr,
        output.loc = "results/MPP_QTL_analysis"
      )
      
    } else{
      
      m <- mpp_proc(
        pop.name = "BCNAM",
        trait.name = traits[i],
        mppData = mppData,
        trait = i,
        Q.eff = models[j],
        threshold = thr,
        output.loc = "results/MPP_QTL_analysis"
      )
      
    }
    
  }
  
}