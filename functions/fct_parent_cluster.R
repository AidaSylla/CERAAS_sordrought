function(mppData, method = NULL, par.clu = NULL,
         w1 = "kernel.exp", w2 = "kernel.unif",
         window, K = 10, simulation.type = "equi",
         simulation.Ng = 50,  simulation.Nrep = 3,
         threshold.quantile = 95, plot = TRUE,
         plot.loc = getwd()){
  
  # check the format of the data
  #################################
  
  if(!is_mppData(mppData)){
    
    stop("'mppData' must be of class ", dQuote("mppData"))
    
  }
  
  # test if correct step in the mppData processing
  
  if(!(mppData$status %in% c('IBD', 'complete'))){
    
    stop("you have to process 'mppData' in a strict order: ",
         "create.mppData, QC.mppData, IBS.mppData, IBD.mppData, ",
         "parent_cluster.mppData. You can only use parent_cluster.mppData ",
         "after create.mppData, QC.mppData, IBS.mppData, and IBD.mppData")
    
  }
  
  # check method
  
  if (is.null(method)){
    
    stop("'method' is not provided")
    
  }
  
  if(!(method %in% c("clusthaplo", "given"))){
    
    stop("'method' must be ", dQuote("clusthaplo"), ' or ', dQuote("given"))
    
  }
  
  if (method == "clusthaplo") {
    
    # Test if clusthaplo is available
    #################################
    
    test <- requireNamespace(package = 'clusthaplo', quietly = TRUE)
    
    if(!test){
      
      stop("the clusthaplo library is not available")
      
    }
    
    # Restore the necessary objects from the mppData object
    ########################################################
    
    haplo.map <- mppData$haplo.map
    consensus.map <- mppData$map[, -3]
    map <- mppData$map
    marker.data <- t(mppData$geno.par.clu)
    parents <- mppData$parents
    
    # cluster the parents
    #####################
    
    # For the step size, determine the value of the largest chromsome
    
    chr.fact <- factor(x = map[, 2], levels = unique(map[, 2]))
    
    step.size <- max(tapply(X = map[, 4], INDEX = chr.fact, FUN = max)) + 100 
    
    p_clu <- parent_cluster(haplo.map = haplo.map, consensus.map = consensus.map,
                            marker.data = marker.data, na.strings = NA,
                            w1 = w1, w2 = w2, step.size = step.size,
                            window = window, K = K,
                            simulation.type = simulation.type,
                            simulation.Ng = simulation.Ng,
                            simulation.Nrep = simulation.Nrep,
                            threshold.quantile = threshold.quantile, plot = plot,
                            plot.loc = plot.loc)
    
    # Check the monomorphic positions
    #################################
    
    par.clu <- parent_clusterCheck(par.clu = p_clu[[1]])
    
    # put the column order as the parents
    
    p_c <- par.clu[[1]]
    p_c <- p_c[, parents]
    
    # Fill the mppData object
    ##########################
    
    mppData$par.clu <- p_c
    
    nb.cl <- apply(X = p_c, MARGIN = 1, FUN = function(x) length(unique(x)))
    
    mppData$n.anc <- mean(nb.cl)
    
    mppData$n.anc.clusthaplo <- p_clu[[2]]
    
    mppData$mono.anc <- par.clu[[2]]
    
    mppData$status <- 'complete'
    
    mppData$geno.off <- NULL
    
    class(mppData) <- c("mppData", "list")
    
    return(mppData)
    
    
  } else { # method = "given"
    
    if(!is.matrix(par.clu)){
      
      stop("'par.clu' argument is not a matrix")
      
    }
    
    if(!is.integer(par.clu)){
      
      stop("'par.clu' is not integer")
      
    }
    
    # list parent
    
    new_par <- colnames(par.clu)
    
    if(!all(new_par %in% mppData$parents)) {
      
      wrong.par <- new_par[!(new_par %in% mppData$parents)]
      pbpar <- paste(wrong.par, collapse = ", ")
      
      message <- sprintf(ngettext(length(wrong.par),
                                  "the following parent %s is not present in 'mppData'",
                                  "the following parents %s are not present in 'mppData'"),
                         pbpar)
      
      stop(message)
      
    }
    
    # list markers
    
    if(!identical(rownames(par.clu), mppData$map[, 1])){
      
      stop("the markers of 'par.clu' and 'mppData' are not identical")
      
    }
    
    # Check monomorphism in par.clu
    ###############################
    
    par.clu <- par.clu[, mppData$parents]
    
    par_clu <- parent_clusterCheck(par.clu = par.clu)
    
    # Calculate the number of ancestral cluster
    ###########################################
    
    nb.cl <- apply(X = par_clu[[1]], MARGIN = 1,
                   FUN = function(x) length(unique(x)))
    
    av.cl <- mean(nb.cl)
    
    mppData$par.clu <- par_clu[[1]]
    
    mppData$n.anc <- av.cl
    
    mppData$mono.anc <- par_clu[[2]]
    
    mppData$status <- 'complete'
    
    mppData$geno.off <- NULL
    
    class(mppData) <- c("mppData", "list")
    
    return(mppData)
    
    
    
  }
  
  
}