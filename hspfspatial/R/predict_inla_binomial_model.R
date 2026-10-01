#' Generate posterior predictions from an INLA binomial spatial model and summarise their distributions.
#'
#' @description
#' Generate posterior predictions from an INLA binomial spatial model and summarise their distributions.
#'
#' @param posterior.samples List of posterior samples returned by INLA.
#'
#' @param mesh INLA mesh object.
#'
#' @param prediction_locations Matrix or data frame of coordinates at which to predict.
#'
#' @param covariates Optional data frame or matrix of covariate values.
#'
#' @return A list containing posterior predictions and their mean, standard deviation, quartiles, and interquartile range.@
#'
#' @export
#'
predict_inla_binomial_model <- function (posterior.samples, mesh, prediction_locations, covariates = NULL) 
{
    nn = length(posterior.samples)
    A.pred <- INLA::inla.spde.make.A(mesh = mesh, loc = prediction_locations)
    mypred <- predict_values(nn, posterior.samples, A.pred = A.pred, 
        covariates = covariates)
    colnames(mypred) = sprintf("posterior_sample_%d", 1:nn)
    pred_mean <- rowMeans(mypred, na.rm = TRUE)
    pred_sd <- apply(mypred, 1, function(x) sd(x, na.rm = TRUE))
    pred_25pct <- apply(mypred, 1, function(x) quantile(x, probs = c(0.25), 
        na.rm = TRUE))
    pred_50pct <- apply(mypred, 1, function(x) quantile(x, probs = c(0.5), 
        na.rm = TRUE))
    pred_75pct <- apply(mypred, 1, function(x) quantile(x, probs = c(0.75), 
        na.rm = TRUE))
    IQR <- pred_75pct - pred_25pct
    return(list(predictions = mypred, mean = pred_mean, sd = pred_sd, 
        q25 = pred_25pct, q50 = pred_50pct, q75 = pred_75pct, 
        iqr = IQR))
}
