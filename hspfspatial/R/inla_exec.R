#' inla_exec
#' @return
#' @export
inla_exec <- function (allModelsList, i) 
{
    formula <- allModelsList[i]
    result <- inla(as.formula(formula), data = inladata, family = "binomial", 
        Ntrials = n, control.predictor = list(A = inla.stack.A(stk), 
            compute = TRUE), control.compute = list(cpo = TRUE, 
            config = TRUE, waic = TRUE, dic = TRUE), list(int.strategy = "eb", 
            diff.logdens = 4), control.fixed = list(prec = myprec, 
            prec.intercept = myprecintercept), verbose = FALSE)
    if (result$ok == FALSE) {
        result <- inla.cpo(result, force = FALSE)
    }
    result_model <- data.frame(Model = as.character(formula), 
        CPO = -1 * mean(log(result$cpo + 0.1), na.rm = TRUE), 
        WAIC = result$waic$waic, DIC = result$dic$dic)
    setTxtProgressBar(mypb, i, title = "Model fit completed", 
        label = i)
    return(result_model)
}
