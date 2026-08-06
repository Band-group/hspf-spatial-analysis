#' makeinlaformula
#' @return
#' @export
makeinlaformula <- function (covariate = NULL) 
{
    if (!is.null(covariate)) {
        myformula0 <- paste("Y ~ -1 + z.intercept + f(z.field, model = spde)")
        myformula <- myformula0
        for (i in 1:length(colnames(covariate))) {
            classcov <- class(covariate[, i])
            covname <- ifelse(classcov == "factor", paste0("f(", 
                colnames(covariate)[i], ", model = \"linear\")"), 
                colnames(covariate)[i])
            myformula <- paste(myformula, covname, sep = " + ")
        }
    }
    else {
        myformula <- paste("Y ~ -1 + z.intercept + f(z.field, model = spde)")
    }
    return(formula(myformula))
}
