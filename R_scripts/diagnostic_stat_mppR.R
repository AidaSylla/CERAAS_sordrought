################################
# build diagnostic stat - mppR #
################################

library(mppR)
library(dplyr)
library(ggplot2)

data("mppData")

# SIM
sim <- mpp_SIM(mppData = mppData)
Q_sim <- QTL_select(sim)
cim <- mpp_CIM(mppData = mppData, cofactors = Q_sim)
plot(cim)
Q_cim <- QTL_select(cim)

QTL_proc <- mpp_proc(mppData = mppData, Q.eff = "par", output.loc = tempdir())
Q_sel <- QTL_proc$QTL

# plot QTL1

# form a data frame object with pheno and ABH marker information

mppData = mppData
QTL = Q_sel
Qpos = 1
trait = "ULA"

QTLeffect_diagnostics <- function(mppData, trait = 1, QTL, Qpos = 1){
  
  # process/check arguments ----
  pheno <- mppData$pheno
  if(is.numeric(trait)){
    trait_id <- colnames(pheno)[trait]
  } else if (is.character(trait)){
    if(!(trait %in% colnames(pheno))){
      message <- paste('trait must be one of:',
                       paste(colnames(pheno), collapse = ", "))
      stop(message)
    } else{
      trait_id <- trait
    }
  }
  
  Q_i <- QTL[Qpos, ]
  d_geno <- mppData$geno.IBD$geno[[Q_i$chr]]
  mk_pos <- which(names(d_geno$map) == Q_i$mk.names)
  
  IBD_proba <- d_geno$prob[, mk_pos, ]
  geno <- c(NA, rep = nrow(IBD_proba))
  for(i in 1:nrow(IBD_proba)){
    geno[i] <- colnames(IBD_proba)[which.max(IBD_proba[i, ])]
  }
  
  # Plot: data process ----
  
  cr_names <- paste(mppData$par.per.cross[, 2], "x", mppData$par.per.cross[, 3])
  names(cr_names) <- mppData$par.per.cross[, 1]
  cross_par_id <- cr_names[mppData$cross.ind]
  
  d_plot <- data.frame(cross = mppData$cross.ind,
                       cross_par_id = cross_par_id,
                       geno = geno,
                       pheno = pheno[, trait_id])
  
  d_plot <- d_plot[!is.na(d_plot$pheno), ]
  
  # Plot: plot ----
  
  p_title <- paste0("Genotype distribution plot - ", trait_id)
  
  p <- ggplot(d_plot, aes(x = geno, y = pheno)) +
  geom_point() +
  geom_smooth(method = "lm", aes(group = 1)) +
  facet_wrap(~ cross_par_id) + ggtitle(p_title)
  
  # summary of N geno per catetories ----
  N_ind_cat <- d_plot |> group_by(cross, cross_par_id, geno) |>
    summarise(N = length(pheno),
              trait_av = mean(pheno), .groups = "drop_last")
  
  # return results ----
  rownames(Q_i) <- NULL
  data = data.frame(Q_i, d_plot)


  return(list(plot = p, data = data, N_ind_cat = N_ind_cat))
  
}


# plot/distribution
Q_diag <- QTLeffect_diagnostics(mppData = mppData, QTL = Q_sel,
                              Qpos = 1,
                              trait = "ULA")

Q_diag$plot

Q_diag$N_ind_cat


# Effect value
par_names <- rownames(QTL_proc$QTL.effects[[1]][[Qpos]])
Qeff <- QTL_proc$QTL.effects[[1]][[Qpos]]$Effect
names(Qeff) <- par_names
Qeff

# % valeur moyenne du trait
trait_av <- mean(mppData$pheno[, trait], na.rm = TRUE)
(Qeff/trait_av) * 100

# R2
R2 <- QTL_proc$R2$part.adj.R2.diff[Qpos]

# Constrasting individual selection ...
data_QTL <- Q_diag$data
