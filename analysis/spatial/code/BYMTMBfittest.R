#TMB version of fitting BYM model

#load packages
library(TMB)
library(Matrix)
library(igraph)
library(spdep)
library(units)
library(sf)


#for testing only###############################################################
mygrid <- grid %>% filter( in_range == 1 );
y_name = sprintf( "%s_+", args$locus );
n_name = sprintf( "%s_N", args$locus );
hbs_columns = grep( "posterior_sample", colnames(hbs), value = T );
model = args$model;
transform = "identity";
link = "logit";
number_of_posterior_samples = args$posterior_samples_per_hbs_sample;
threads = args$threads;
################################################################################

TMBfitbym_to_posterior_samples <- function(
    mygrid, hbs, pf,
    y_name = "Pfsa1_+",
    n_name = "Pfsa1_N",
    hbs_columns = "posterior_mean",
    model = "bym2", # or "iid" or "norandom" or "besag"
    transform = "identity",
    link = "logit",
    number_of_posterior_samples = 100,
    threads = 1
) {
  countrydfi = (
    mygrid
    %>% dplyr::left_join( pf, by = "polygon_id" )
  )
  print( countrydfi )
  countrydfi$Y = countrydfi[[y_name]]
  countrydfi$n = countrydfi[[n_name]]
  
  ### remove polygon if missing response or sample size
  countrydfi <- countrydfi %>% dplyr::filter(!is.na(Y) & !is.na(n))
  countrydfi <- sf::st_make_valid(countrydfi)
  echo( "++ data for fitting is:" )
  print( countrydfi, n = 50 )
  
  ############################################################
  hbs = hbs[ match( countrydfi$polygon_id, hbs$polygon_id ), ]
  print( dim( countrydfi ))
  print( dim( hbs ))
  
  #check if redundant polygon_id?
  stopifnot(
    length(
      countrydfi
      %>% dplyr::group_by(polygon_id)
      %>% dplyr::filter(n() > 1)
      %>% dplyr::pull(polygon_id)
      %>% unique()
    ) == 0
  )
  
  #RINLA needs ID from 1 to ...otherwise leads to issue during fitting process
  countrydfi$ID <- 1:nrow(countrydfi)
  
  transform.fn = get( transform )
  #formula for BYM model with pc priors (without  F.E., to be updated)
  if( model == 'norandom' ) {
    myformula <- (
      Y ~ -1
      + intercept
      + transform.fn(HbAS_or_SS)
    )
    
  } else if( model == 'iid' ) {
    myformula <- (
      Y ~ -1
      + intercept
      + transform.fn(HbAS_or_SS)
      + f( ID, model = "iid" )
    )
    
  } else if( model == 'bym2' ) {
    myformula <- (
      Y ~ -1
      + intercept
      + transform.fn(HbAS_or_SS)
      + f( ID, model = model, graph = g, hyper = prior, scale.model = TRUE, constr = TRUE )
    )
 
  } else if( model == 'besag' ) {
    myformula <- (
      Y ~ -1
      + intercept
      + transform.fn(HbAS_or_SS)
      + f( ID, model = model, graph = g )
    )
    } else {
    stop( sprintf( "Unrecognised model \"%s\".  (I only support 'norandom', 'besag', 'bym2' or 'iid' currently.)", model ))
    }
  
  #C++ file to be compiled ################
  compile(paste0("code/",model,".cpp"))
  dyn.load(dynlib(paste0("code/",model)))
  #dyn.unload(dynlib(paste0("code/",model)))
  #########################################

  fitted.parameters = tibble()
  sampled.parameters = tibble()
  summary = tibble()
  
  n = nrow(countrydfi)
  
  #if spatial term in the model
  #Define spatial matrix and all the necessary for running TMB BYM
  #update this if necessary
  if(model %in% c('besag','bym2')) {
  nb <- spdep::poly2nb(countrydfi,queen = TRUE)#,snap=mysnap)
  # Find nodes without neighbors
  mstconnect <- function(polys, nb, distance="centroid"){
    if(distance == "centroid"){
      coords = sf::st_coordinates(sf::st_centroid(sf::st_geometry(polys)))
      dmat = as.matrix(dist(coords))
    } else if(distance == "polygon"){
      dmat = sf::st_distance(polys) + units::set_units(1000, "m") # offset for adjacencies
      diag(dmat) = 0 # no self-intersections
    }else{
      stop("Unknown distance method")
    }
    
    gfull = igraph::graph_from_adjacency_matrix(dmat, weighted=TRUE, mode="undirected")
    gmst = igraph::mst(gfull)
    #gmst = gfull
    edgemat = as.matrix(igraph::as_adj(gmst))
    edgelistw = spdep::mat2listw(edgemat,style="M")
    edgenb = edgelistw$neighbour
    attr(edgenb,"region.id") = attr(nb, "region.id")
    allnb = spdep::union.nb(nb, edgenb)
    return(allnb)
  }
  #slow if polygon distance is used#####################
  nball <- mstconnect(countrydfi, nb, distance="centroid")
  ######################################################
  # check the network
  #plot(st_geometry(countrydfi), border = "grey");
  #coords = sf::st_coordinates(sf::st_centroid(sf::st_geometry(countrydfi)))
  #plot(nball,coords, add = T,lwd=1, col="red")
  #plot(nb, coords, add = T, col="blue",lwd=0.5)
  ######################################################
  
  td = tempdir()
  tempfile = sprintf( "%s/%s", td, "countrydfi.adj" )
  spdep::nb2INLA( tempfile, nball)
  g <- INLA::inla.read.graph(filename = tempfile )
  adj_matrix <- INLA::inla.graph2matrix(g)
  } 
  #scale Q matrix
  # Q = -inla.graph2matrix(g)
  # diag(Q) = 0
  # diag(Q) = -rowSums(Q)
  # n = dim(Q)[1]
  # Q.scaled <- inla.scale.model(Q,constr = list(A = matrix(1, 1, n), e=0))
 
  #TMB BYM parameters (initial values)
  TMBpara <- list(
    intercept = 0.1, 
    HbAS_or_SS = 1, 
    u = rep(0.1, n),  # Assuming 'y' is your response vector
    v = rep(0.1, n-1),  # v has dimension n-1 because it is forced to sum to 0 in the TMB model
    log_tau_u = 0.1,  # Neutral start
    log_tau_v = 1.0   # Neutral start
  )
  
  #loop to run regression for each posterior sample of HbS
  for( sample in hbs_columns ) {
   #sample <- hbs_columns[1] #for test
    regression.data <- data.frame(
      countrydfi,
      intercept = rep(1, length(countrydfi$Y)),
      HbAS_or_SS = hbs[[sample]]^2 + 2*hbs[[sample]]*(1-hbs[[sample]])
    )
   
    if(model %in% c('bym2')) {
    data <- list(y = regression.data$Y, N = regression.data$n, x = regression.data$HbAS_or_SS, 
                 adj_matrix = adj_matrix)
    myrandom <- c('u','v')
    }
    if(model %in% c('besag')) {
      data <- list(y = regression.data$Y, N = regression.data$n, x = regression.data$HbAS_or_SS, 
                   adj_matrix = adj_matrix)
      myrandom <- c('v')
    }
    if(model %in% c('iid')) {
      data <- list(y = regression.data$Y, N = regression.data$n, x = regression.data$HbAS_or_SS, 
                   adj_matrix = adj_matrix)
      myrandom <- c('u') 
    }
      if(model %in% c('norandom')) {
        data <- list(y = regression.data$Y, N = regression.data$n, x = regression.data$HbAS_or_SS, 
                     adj_matrix = adj_matrix)
        myrandom <- NULL
      }
    # Build and optimize model
    obj <- TMB::MakeADFun(
      data = data,
      parameters = TMBpara,
     random = myrandom, #u: iid term, v: spatial term, considered 'random effects'
      DLL = model,
     inner.control = list(
       maxit = 10000,           # Increase maximum iterations
       tol = 1e-8,             # Tolerance for convergence
       trace = TRUE,           # Print progress
       step.tol = 1e-12,       # Step tolerance
       mgcmax = 1e+20,         # Maximum gradient component
       sir = TRUE,             # Use saddle point approximation if needed
       newton = TRUE          # Avoid Newton method if causing issues
     )
    )
    #various options to fit the model (optimization procedures)
    #Option 1: using nlminb
    opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1000, iter.max = 1000))
    
    #Option 2: using optim (different variation: BFGS performs better when dimension of parameter is high.
    #opt <- optim(obj$par, obj$fn, obj$gr, method = "BFGS",control=list(maxit = 25000,ndeps=0.001))     # 
    
    # Get estimates and uncertainties
    # This will need to be updated (using MCMC, other? to get uncertainty)
    report <- sdreport(obj)
    #summary(report, "fixed") 
    #summary(report, "random") 
   
    #INLA Version###############################################################
    # fit <- INLA::inla(
    #   myformula,
    #   family = "binomial",
    #   control.family = list( control.link = list( model = link )),
    #   data = regression.data,
    #   Ntrials = n, # this is specific to binomial as we need to tell it the number of examined
    #   control.predictor = list(compute = TRUE), # compute gives you the marginals of the linear predictor
    #   control.compute = list(return.marginals.predictor=TRUE, waic = TRUE, cpo = TRUE, mlik=TRUE, config = TRUE), # model diagnostics and config = TRUE gives you the GMRF,mlik = TRUE to compute marg.likelihood
    #   control.inla = list(strategy = "laplace", npoints = 21),#better approximation and increase evaluation points
    #   #list(int.strategy = "grid", diff.logdens = 4),#to improve CPO computation
    #   verbose = FALSE,
    #   num.thread = threads
    # )
    # #summary of results
    # s = summary(fit)
    # print( s )
    # # We store BOTH the parameter fits
    # and also a sample of posterior parameters from the model
    # for later visualisation
    # fitted.parameters = bind_rows(
    #   fitted.parameters,
    #   bind_cols(
    #     hbs.sample = sample,
    #     model = model,
    #     parameter = c( rownames(s$fixed), rownames(s$hyperpar) ),
    #     rbind( s$fixed[,1:6], s$hyperpar[,1:6] )
    #   )
    # )
    # summary = bind_rows(
    #   summary,
    #   tibble(
    #     hbs.sample = sample,
    #     model = model,
    #     cpo = -1*mean( log(fit$cpo$cpo+0.1), na.rm = TRUE),
    #     waic = fit$waic$waic,
    #     marginal_ll_integration = s$mlik[1],
    #     marginal_ll_gaussian = s$mlik[2]
    #   )
    # )
    
    #sample from posterior distribution of INLA    
    # posterior.parameters = inla.posterior.sample( number_of_posterior_samples, fit )
    # sampled.parameters = bind_rows(
    #   sampled.parameters,
    #   bind_cols(
    #     hbs.sample = sample,
    #     model = model,
    #     intercept = sapply(
    #       posterior.parameters,
    #       function(x) {
    #         x$latent[grep('intercept', rownames(x$latent)),1]
    #       }
    #     ),
    #     beta = sapply(
    #       posterior.parameters,
    #       function(x) {
    #         x$latent[grep('HbAS_or_SS', rownames(x$latent)),1]
    #       }
    #     )
    #   )
    # )
    ############################################################################
    #TMB version################################################################
    #create approx. 95 CI and mode
    TMBfixfit <- data.frame(summary(report, "fixed"))
    colnames(TMBfixfit) <- c('mean','sd')
    TMBfixfit$`0.025quant` <- TMBfixfit$mean - 1.96 * TMBfixfit$sd
    TMBfixfit$`0.5quant` <- TMBfixfit$mean
    TMBfixfit$`0.975quant` <-  TMBfixfit$mean + 1.96 * TMBfixfit$sd
    TMBfixfit$mode <- TMBfixfit$mean
        
    fitted.parameters = bind_rows(
      fitted.parameters,
      bind_cols(
        hbs.sample = sample,
        model = model,
        parameter = c( names(report$par.fixed)),
        TMBfixfit)
      )
    summary = bind_rows(
      summary,
      tibble(
        hbs.sample = sample,
        model = model,
        cpo = NA,
        waic = report$gradient.fixed[1],#not waic but max gradient component
        marginal_ll_integration = NA,
        marginal_ll_gaussian = NA
      )
    )
    #replace with TMB distribution sampling#####################################
    ############################################################################
    #posterior.parameters = inla.posterior.sample( number_of_posterior_samples, fit )
    #a very rough sampling (before implementing TMBstan or other)
    hbs.sample <- rnorm(n=number_of_posterior_samples, 
          mean=fitted.parameters[fitted.parameters$parameter=='HbAS_or_SS',]$mean,
          sd=fitted.parameters[fitted.parameters$parameter=='HbAS_or_SS',]$sd)
    intercept.sample <- rnorm(n=number_of_posterior_samples, 
                        mean=fitted.parameters[fitted.parameters$parameter=='intercept',]$mean,
                        sd=fitted.parameters[fitted.parameters$parameter=='intercept',]$sd)
    
    sampled.parameters <- data.frame(
      hbs.sample = rep(sample,number_of_posterior_samples),
      model = rep(model,number_of_posterior_samples),
      intercept = intercept.sample,
      beta = hbs.sample)
   
    sampled.parameters <- as_tibble(sampled.parameters)
    
    ############################################################################
   
    echo( "... ++ Ok, successfully fit model for %s..\n", sample )
  }
  # fix parameter name for the transform
  fitted.parameters$parameter = gsub( "^transform[.]fn[(]", sprintf( "%s(", transform ), fitted.parameters$parameter )
  return(
    list(
      model = model,
      data = countrydfi,
      transform = transform,
      link = link,
      prior = prior,
      allele = y_name,
      fitted.parameters = fitted.parameters,
      sampled.parameters = sampled.parameters,
      summary = summary
    )
  )
}
