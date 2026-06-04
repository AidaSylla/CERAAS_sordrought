library(mppR)
#> mppR has been developed with the support of
#> Wageningen University
#> KWS SAAT SE & Co. KGaA
#> Swiss National Science Foundation
#> (Grant Postdoc.Mobility P500PB_203030)

#Data
data(mppData_GE)
design_connectivity(par_per_cross = mppData_GE$par.per.cross)

#Simple interval mapping (SIM)
SIM <- mppGE_SIM(mppData = mppData_GE, trait = c('DMY_CIAM', 'DMY_TUM'), ref_par = 'UH007')
plot(x = SIM)

#Cofactors selection
cofactors <- QTL_select(Qprof = SIM, threshold = 4, window = 50)

#Composite interval mapping (CIM)
CIM <- mppGE_CIM(mppData = mppData_GE, trait = c('DMY_CIAM', 'DMY_TUM'),
                 cofactors = cofactors, window = 20, VCOV_data = "unique")

##plot(x = CIM)
#QTLs selection
QTL <- QTL_select(Qprof = CIM, threshold = 4, window = 20)

#QTLs effect and significance estimations
Qeff <- QTL_effect_GE(mppData = mppData_GE, trait = c('DMY_CIAM', 'DMY_TUM'),
                      QTL = QTL)
Qeff$QTL_1
#>                 Effect Std.dev  Wald df        p.val sign
#> QTL_1_UH007_E1      NA      NA    NA NA           NA <NA>
#> QTL_1_D152_E1    0.979   3.107  0.10  1 0.7526802395     
#> QTL_1_F03802_E1 -4.632   2.265  4.18  1 0.0408643496    *
#> QTL_1_F2_E1     -3.170   3.523  0.81  1 0.3683083123     
#> QTL_1_F283_E1   -7.691   2.189 12.35  1 0.0004412536  ***
#> QTL_1_DK105_E1  -9.861   2.615 14.22  1 0.0001626113  ***
#> QTL_1_UH007_E2      NA      NA    NA NA           NA <NA>
#> QTL_1_D152_E2    2.351   3.010  0.61  1 0.4346726625     
#> QTL_1_F03802_E2  1.921   2.173  0.78  1 0.3764976065     
#> QTL_1_F2_E2     -3.200   3.386  0.89  1 0.3445433645     
#> QTL_1_F283_E2   -0.561   2.110  0.07  1 0.7904495664     
#> QTL_1_DK105_E2  -1.559   2.513  0.38  1 0.5350108782

#Ici on utilise F2 comme parent réccurrent
Qeff <- QTL_effect_GE(mppData = mppData_GE, trait = c('DMY_CIAM', 'DMY_TUM'),
                      QTL = QTL, ref_par = 'F2')
Qeff$QTL_1
#>                 Effect Std.dev Wald df     p.val sign
#> QTL_1_UH007_E1   3.170   3.523 0.81  1 0.3683083     
#> QTL_1_D152_E1    4.149   4.697 0.78  1 0.3771534     
#> QTL_1_F03802_E1 -1.462   4.188 0.12  1 0.7269845     
#> QTL_1_F2_E1         NA      NA   NA NA        NA <NA>
#> QTL_1_F283_E1   -4.522   4.148 1.19  1 0.2756179     
#> QTL_1_DK105_E1  -6.691   4.387 2.33  1 0.1272424     
#> QTL_1_UH007_E2   3.200   3.386 0.89  1 0.3445434     
#> QTL_1_D152_E2    5.552   4.530 1.50  1 0.2203929     
#> QTL_1_F03802_E2  5.122   4.023 1.62  1 0.2029702     
#> QTL_1_F2_E2         NA      NA   NA NA        NA <NA>
#> QTL_1_F283_E2    2.640   3.989 0.44  1 0.5081781     
#> QTL_1_DK105_E2   1.641   4.217 0.15  1 0.6970780

#QTLs R2 calculation
QR2 <- QTL_R2_GE(mppData = mppData_GE, trait = c('DMY_CIAM', 'DMY_TUM'), QTL = QTL)

#la contribution des QTL
QR2$glb.adj.R2
#> [1] 6.328735
# la contribution de chaque QTL
QR2$part.adj.R2.diff
#>       Q1       Q2 
#> 2.563632 3.582575


#QTL profile plot
plot(x = CIM)


#Whole-genome genetic effect significance plot
plot_allele_eff_GE(mppData = mppData_GE, nEnv = 2,
                   EnvNames = c('CIAM', 'TUM'), Qprof = CIM,
                   QTL = QTL, text.size = 14)

#mppGE_proc: wrapper function
MPP_GE_QTL <- mppGE_proc(pop.name = 'EUNAM', trait.name = 'DMY',
                         mppData = mppData_GE, trait = c('DMY_CIAM', 'DMY_TUM'),
                         n.cores = 1, verbose = FALSE, output.loc = tempdir())

#QTL by environment effect determination
Qeff <- QTL_effect_main_QEI(mppData = mppData_GE,
                            trait = c('DMY_CIAM', 'DMY_TUM', 'DMY_INRA_P', 'DMY_KWS'),
                            env_id = c('CIAM', 'TUM', 'INRA', 'KWS'),
                            QTL = QTL)

Qeff$Q_sign$QTL2
#>      par  logP_main  logP_QxE
#> 1  UH007         NA        NA
#> 2     F2 1.05793090 3.1420765
#> 3   D152 0.96910440 0.2871339
#> 4  DK105 0.05370229 0.2525703
#> 5   F283 0.36207766 3.9077148
#> 6 F03802 0.01518149 4.6439537

Qeff$Q_eff$QTL2

#QTLxEC effect estimation
# provide the environmental covariate information as a matrix
EC <- matrix(c(180, 310, 240, 280), 4, 1)
rownames(EC) <- c('CIAM', 'TUM', 'INRA', 'KWS')
colnames(EC) <- 'cum_rain'

Qeff <- QTL_effect_main_QxEC(mppData = mppData_GE,
                             trait = c('DMY_CIAM', 'DMY_TUM', 'DMY_INRA_P', 'DMY_KWS'),
                             env_id = c('CIAM', 'TUM', 'INRA', 'KWS'),
                             QTL = QTL, EC = EC, Qmain_QEI = Qeff, thre_QTL = 1.301)

Qeff$Qeff_EC$QTL2

#Plot QTL x EC sensitivity slopes
plot_QxEC(Qeff, EC = EC, env_id = c('CIAM', 'TUM', 'INRA', 'KWS'), 
          QTL = 2, EC_id = 'cum rain', trait_id = 'DMY')
