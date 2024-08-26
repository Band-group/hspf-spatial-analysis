# BYM
#setwd("D:/OneDrive/MOCHIALL/MOCHI/PROJECT/MED/MED2_HBSPF/hspf-spatial-analysis/analysis/spatial/")
#getwd()

# Note on generalising this:
# - input grid type (hexagon/squares)
# - properly use the posterior sample of HbS (only one sample used so far)
# - Pf only 1st allele used
# - No covariates / fixed effects at the mo
# - 

library(spdep)
library(INLA)
library(dplyr)

#get polygon data 
grid <- readRDS('output/grids/grid-type=hexagon-size=1-division=none.rds')
#get hbs data
hbs <- readr::read_tsv('output/HbSsensitivity/fixed-r0=10.0-sigma0=1.0-fc=none/aggregated/grid-type=hexagon-size=1-division=none.tsv')
#get Pf data
# Some hexagons have data from two source studies, aggregate as we load
pf = (
  readr::read_tsv('output/HbSsensitivity/pf/aggregated/grid-type=hexagon-size=1-division=none.tsv')
  %>% group_by( polygon_id )
  %>% summarise(
    `Pfsa1_+` = sum(`Pfsa1_+`),
    Pfsa1_N = sum( Pfsa1_N ),
    `Pfsa2_+` = sum( `Pfsa2_+` ),
    Pfsa2_N = sum( Pfsa2_N ),
    `Pfsa3_+` = sum(`Pfsa3_+`),
    Pfsa3_N = sum( Pfsa3_N ),
    `Pfsa4_+` = sum( `Pfsa4_+` ),
    Pfsa4_N = sum( Pfsa4_N )
  )
)

#extract hbs
#do it for one sample only (need to be done for all samples...)
kepvar <- c('polygon_id','posterior_sample_1')
myhbs <- myhbs[kepvar]
colnames(myhbs) <- c('polygon_id','HbS')
#add HbS data
countrydfi <- (
  countrydfi
  %>% dplyr::left_join(myhbs, by = c("polygon_id"))
  %>% dplyr::left_join(pf, by = "polygon_id")
)

######TO be done?
#remove polygon if missing covariates
#seems ok without it
#countrydfi <- countrydfi %>%
#  dplyr::filter(!is.na(HbS))
stopifnot( length(which(is.na( countrydfi$HbS ))) == 0 )

#extract pf (only for pfsa1 need to be done for other alleles)

priors = tibble::tribble(
  ~item, ~threshold, ~p_greater_than,
  "pc.prec", 0.5 / 0.31, 0.01,
  "pc", 0.5, 2/3
)

fitbym <- function(
  grid,
  hbs,
  pf,
  locus = 1,
  model = "bym2", # or "iid" or "norandom" or "besag"
  # TODO: Priors currently only work with bym2 model, fix this.
  prior = list(
    prec = list(
      prior = "pc.prec",
      param = c(0.5 / 0.31, 0.01)),
    phi = list(
      prior = "pc",
      param = c(0.5, 2 / 3))
  ),
  number_of_posterior_samples = 25
) {
  countrydfi = (
    grid
    %>% dplyr::left_join( pf, by = "polygon_id" )
  )
  pfsa_plus_name = sprintf( 'Pfsa%d_+', locus )
  pfsa_N_name = sprintf( 'Pfsa%d_N', locus )
  countrydfi$Y = countrydfi[[pfsa_plus_name]]
  countrydfi$n = countrydfi[[pfsa_N_name]]

  hbs = hbs[ match( countrydfi$polygon_id, hbs$polygon_id ), ]

  ###remove polygon if missing response or sample size
  countrydfi <- countrydfi %>% dplyr::filter(!is.na(Y) & !is.na(n))
  ############################################################

  #check if redundant polygon_id?
  stopifnot(
    length(
      countrydfi %>%
        dplyr::group_by(polygon_id) %>%
        dplyr::filter(n() > 1) %>%
        dplyr::pull(polygon_id) %>%
        unique()
    ) == 0
  )

  #RINLA needs ID from 1 to ...otherwise leads to issue during fitting process
  countrydfi$ID <- 1:nrow(countrydfi)

  #set pc prior for spatial and unstructured term
  # prior should have prec and phi entries
  # each with prior and param entries
  #prior <- list(
  #  prec = list(
  #    prior = "pc.prec",
  #    param = c(0.5 / 0.31, 0.01)),
  #  phi = list(
  #    prior = "pc",
  #    param = c(0.5, 2 / 3))
  #)

  #create adjacent matrix
  nb <- spdep::poly2nb(countrydfi)
  td = tempdir()
  tempfile = sprintf( "%s/%s", td, "countrydfi.adj" )
  spdep::nb2INLA( tempfile, nb)
  g <- INLA::inla.read.graph(filename = tempfile )

  #formula for BYM model with pc priors (without  F.E., to be updated)
  if( model == 'norandom' ) {
    myformula <- (
      Y ~ -1
      + y.intercept
      + HbAS_or_SS
    )
  } else if( model == 'iid' ) {
    myformula <- (
      Y ~ -1
      + y.intercept
      + HbAS_or_SS
      + f( ID, model = "iid" )
    )
  } else if( model == 'bym2' ) {
    myformula <- (
      Y ~ -1
      + y.intercept
      + HbAS_or_SS
      + f( ID, model = model, graph = g, hyper = prior, scale.model = TRUE, constr = TRUE )
    )
  } else {
    myformula <- (
      Y ~ -1
      + y.intercept
      + HbAS_or_SS
      + f( ID, model = model, graph = g )
    )
  }

  fitted.parameters = tibble()
  sampled.parameters = tibble()
  # Pray this works
  for( sample in grep( "posterior_sample", colnames(hbs))) {
    regression.data <- data.frame(
      countrydfi,
      y.intercept = rep(1, length(countrydfi$Y)),
      HbAS_or_SS = hbs[[sample]]^2 + 2*hbs[[sample]]*(1-hbs[[sample]])
    )
    #update this
    fit <- INLA::inla(
      myformula,
      family = "binomial",
      data = regression.data,
      Ntrials = n, # this is specific to binomial as we need to tell it the number of examined
      control.predictor = list(compute = TRUE), # compute gives you the marginals of the linear predictor
      control.compute = list(return.marginals.predictor=TRUE, waic = TRUE, cpo = TRUE, config = TRUE), # model diagnostics and config = TRUE gives you the GMRF
      control.inla = list(strategy = "laplace", npoints = 21),#better approximation and increase evaluation points
      #list(int.strategy = "grid", diff.logdens = 4),#to improve CPO computation
      verbose = FALSE,
      num.thread = 1
    )
  #summary of results
  s = summary(fit)
  print( s )
  fitted.parameters = bind_rows(
    fitted.parameters,
    bind_cols(
      hbs.sample = sample,
      model = model,
      parameter = rownames(s$fixed),
      s$fixed
    )
  )
  posterior.parameters = inla.posterior.sample( number_of_posterior_samples, fit )
  sampled.parameters = bind_rows(
    sampled.parameters,
    bind_cols(
      hbs.sample = sample,
      model = model,
        intercept = sapply(
          posterior.parameters,
          function(x) {
            x$latent[grep('y.intercept', rownames(x$latent)),1]
          }
        ),
        beta = sapply(
          posterior.parameters,
          function(x) {
            x$latent[grep('HbAS_or_SS', rownames(x$latent)),1]
          }
        )
    )
  )
  }
  return(
    list(
      model = model,
      prior = prior,
      fitted.parameters = fitted.parameters,
      sampled.parameters = sampled.parameters
    )
  )
}

################################################################################
#Some maps for double check
# countrydfi$RR <- bym$summary.fitted.values[, "mean"]
# countrydfi$LL <- bym$summary.fitted.values[, "0.025quant"]
# countrydfi$UL <- bym$summary.fitted.values[, "0.975quant"]
# mapsf <- sf::st_as_sf(countrydfi)
# library(ggplot2)
# gRR <- ggplot(mapsf) + geom_sf(aes(fill = RR)) +
#   scale_fill_gradient2(
#     midpoint = 1, low = "blue", mid = "white", high = "red",
#     limits = c(0.7, 1.5)
#   ) +
#   theme_bw()
# gLL <- ggplot(mapsf) + geom_sf(aes(fill = LL)) +
#   scale_fill_gradient2(
#     midpoint = 1, low = "blue", mid = "white", high = "red",
#     limits = c(0.7, 1.5)
#   ) +
#   theme_bw()
# gUL <- ggplot(mapsf) + geom_sf(aes(fill = UL)) +
#   scale_fill_gradient2(
#     midpoint = 1, low = "blue", mid = "white", high = "red",
#     limits = c(0.7, 1.5)
#   ) +
#   theme_bw()
# 
# library(cowplot)
# plot_grid(gRR, gLL, gUL, ncol = 1)
