#' predict_values
#' @return
#' @export
predict_values <- function (nn, posterior.samples, A.pred, covariates = NULL, link.function = stats::plogis) 
{
    pred <- matrix(NA, nrow = dim(A.pred)[1], ncol = nn)
    for (i in 1:nn) {
        field <- posterior.samples[[i]]$latent[grep("z.field", 
            rownames(posterior.samples[[i]]$latent)), ]
        intercept <- posterior.samples[[i]]$latent[grep("z.intercept", 
            rownames(posterior.samples[[i]]$latent)), ]
        if (is.null(covariates)) {
            lp <- drop(A.pred %*% field) + intercept
        }
        else {
            beta <- NULL
            linpred <- list()
            k <- ncol(covariates)
            for (j in 1:k) {
                beta[j] <- posterior.samples[[i]]$latent[grep(names(covariates)[j], 
                  rownames(posterior.samples[[i]]$latent)), ]
                linpred[[j]] <- beta[j] * covariates[, j]
            }
            linpred <- Reduce("+", linpred)
            lp <- drop(A.pred %*% field) + intercept + linpred
        }
        pred[, i] <- link.function(as.numeric(lp))
    }
    return(pred)
}
