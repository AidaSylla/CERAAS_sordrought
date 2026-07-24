#####################################################
# IBD.mppData with direct insertion of the ABH data #
#####################################################

# mppData
# type = 'BCsFt'
# F.gen = 3
# BC.gen = 1
# type.mating = "selfing"
# error.prob = 1e-04
# map.function = "kosambi"

IBD.mppData_BCNAM <- function(mppData, type = 'BCsFt', F.gen = 3,
                              BC.gen = 1, type.mating = "selfing",
                              error.prob = 1e-04, map.function = "kosambi"){
  
  ##### necessary objects from the mppData object ######
  
  geno.off <- mppData$geno.off
  geno.par <- mppData$geno.par
  pheno <- mppData$pheno
  map <- mppData$map
  
  #### type of population ####
  
  if (type == "F") {
    
    type.pop <- paste0("F", "(n = ", F.gen,")")
    
  } else if (type == "BC") {
    
    type.pop <- paste0("Back-cross ", "(n = ", BC.gen,")")
    
  } else if (type == "DH") {
    
    type.pop <- "Double haploid"
    
  } else if (type == "RIL") {
    
    type.pop <- "Recombinant inbred line"
    
    if(type.mating == "selfing") {
      
      type.pop <- paste (type.pop, "by selfing")
      
    } else {
      
      type.pop <- paste (type.pop, "by sibling mating")
      
    }
    
  } else if (type == 'BCsFt'){
    
    type.pop <- paste0('Back-cross followed by selfing ', '(',
                       paste0('BC', BC.gen, 'F', F.gen), ')')
    
  }
  
  #### number of allele class ####
  
  if ((type == "BC") | (type == "DH") | (type == "RIL")) {
    
    n.zigo <- 2
    
  } else if ((type == "F")| (type == "BCsFt")) {
    
    n.zigo <- 3
    
  }
  
  
  #### Form the cross object  ####
  
  # convert heterozygous into missing for RIL and DH population
  
  if(type %in% c('RIL', 'DH')){
    
    geno.off[geno.off == 'H'] <- NA
    
  }
  
  # format data to form a cross object
  
  chr.info <- t(map[, 2:3])
  colnames(chr.info) <- map[, 1]
  
  geno.aug <- rbind(chr.info, geno.off)
  
  n_pheno <- dim(pheno)[2]
  
  empty_mat <- matrix("", nrow = 2, ncol = n_pheno)
  
  trait.aug <- rbind(empty_mat, pheno)
  geno.aug <- cbind(trait.aug, geno.aug)
  
  # Export the data in a .csv file in a temporary file.
  
  tmp <- tempfile(fileext = "Cross_object.csv")
  
  write.csv(geno.aug, file = tmp, row.names = FALSE)
  
  rm(geno.aug)
  
  # form a R/qtl cross object reading the data using the specifiec type of
  # population
  
  if (type == "F") {
    
    cross.object <- read.cross("csv", , tmp, F.gen = F.gen)
    
    
  } else if (type == "BC") {
    
    cross.object <- read.cross("csv", , tmp, BC.gen = BC.gen)
    
    
  } else if (type == "RIL") {
    
    # need to read the object as a backcross
    
    cross.object <- read.cross("csv", file = tmp, genotypes = c("A", "B"),
                               alleles = c("A", "B"))
    
    # then convert it following the type of mating
    
    if (type.mating == "selfing") {
      
      cross.object <- convert2riself(cross.object)
      
    }
    
    if (type.mating == "sib.mat") {
      
      cross.object <- convert2risib(cross.object)
      
    }
    
  } else if (type == "DH") {
    
    # need to read the object as a backcross
    
    cross.object <- read.cross("csv", file = tmp, genotypes = c("A", "B"),
                               alleles = c("A", "B"))
    
    class(cross.object)[1] <- "dh"
    
  } else if (type == "BCsFt"){
    
    cross.object <- read.cross("csv", , tmp, F.gen = F.gen, BC.gen = BC.gen)
    
  }
  
  # 6. Compute the IBD probabilities
  ##################################
  
  # make sure no need to provide stepwidth
  
  cross.object <- calc.genoprob(cross.object, step = 0,
                                error.prob = error.prob,
                                map.function = map.function)
  
  # delete the temporary directory
  
  unlink(tmp)
  
  # 7. transform the map and geno.par
  ###################################
  
  chr.ind <- factor(x = map[, 2], levels = unique(map[, 2]))
  
  new.map <- data.frame(map[, 1:2], sequence(table(chr.ind)),
                        map[, 3], stringsAsFactors = FALSE)
  
  colnames(new.map) <- c("mk.names", "chr", "pos.ind", "pos.cM")
  
  ### geno.par
  
  geno.par.new <- t(geno.par)
  
  geno.par.new <- data.frame(new.map, geno.par.new, stringsAsFactors = FALSE,
                             check.names = FALSE)
  
  
  # 8. fill the mppData object
  #############################
  
  mppData$geno.IBD <- cross.object
  
  mppData$type <- type.pop
  
  mppData$n.zigo <- n.zigo
  
  mppData$map <- new.map
  
  mppData$geno.par <- geno.par.new
  
  mppData$status <- 'IBD' 
  
  class(mppData) <- c("mppData", "list")
  
  return(mppData)
  
}