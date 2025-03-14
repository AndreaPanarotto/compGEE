#' GEE fit of compositional time series
#'
#' \code{compGeeFit} returns the parameter estimates for a GEE model, based on the compositional logit model by Firth and Sammut (2023), for longitudinal compositional data analysis. The time series are modeled directly on the simplex, without the requirement of logratio transformations. The model takes into account information about the total, so it is not necessary to provide normalize data.
#'
#' @param Y_list A list of N time series. Each element of the list should be a matrix where each row represents a time observation of dimension D.
#' @param X_list A list of covariates relative to the time series. Each element of the list should be a matrix where each row represents the covariates (of dimension P) for a time observation. See `Details`.
#' @param Beta_0 Starting point for the regression parameters. Should be a matrix of size D x P.
#' @param rho Initial value for the autoregressive correlation.
#' @param Phi Initial value for the error matrix.
#' @param jack Boolean; whether to compute the 1-step jackknife estimate and variance of the estimate from Lipsitz et al. (1990).
#' @param criterion Identifiability criterion for the regression parameters. See `Details`.
#' @param max.iter Maximum number of iterations of the estimation procedure.
#' @param eps Tolerance in update difference for convergence.
#' @param burnin Number of iterations after which the update length of the Fisher's scoring step starts decaying. Default is at \code{max.iter} (= no decay).
#' @param decay_factor Parameter to regulate the strength of the decay.
#' @param verbose A non-negative integer to indicate if and how frequently print the current parameters as the iterations proceed.
#' @details
#' If \code{Y_list} is of length N, so should be \code{X_list.} Then, to each matrix of size T x D in \code{Y_list} should correspond a matrix of size T x P in \code{X_list}: each row represents the covariates for the corresponding row in the observations in \code{Y_list}.
#'
#' Each element \eqn{\beta_{dp}} represents the effect of p-th covariate on the d-th part. For identifiability reasons, we need to impose some conditions on \eqn{\boldsymbol{\beta}}. The reference constraint simply puts the last row of \eqn{\boldsymbol{\beta}}, that is, \eqn{\beta_{Dp} = 0} for all p. In practice, the last component is seen as the reference to which the effect on the others should be compared. The sum constraint, instead, assumes \eqn{\sum_d \beta_{dp} = 0} for all p, so we are imposing the average effect of each covariate to be 0 across the parts.
#'
#' The stopping criterion is on the maximum absolute difference between the elements consecutive instances of \eqn{\boldsymbol{\beta}}.
#'
#' If \code{verbose} is 0, nothing is printed during the computation. If it is a positive integer, the current parameters, the step size, and the difference in the regression parameters from the previous iteration is printed every \code{verbose} iterations.
#'
#' @return A list containing:
#' \describe{
#' \item{Beta_hat}{Matrix of final regression parameters.}
#' \item{rho_hat}{Estimated autocorrelation.}
#' \item{Phi_hat}{Estimated variance matrix of the errors}
#' \item{conv}{A boolean telling whether the algorithm reached convergence.}
#' \item{diff}{Final difference between updates of the regression parameters.}
#' \item{out_diff}{List of differences between updates.}
#' \item{n_iter}{Number of performed iterations.}
#' \item{time}{Computational time.}
#' \item{alg.params}{Some of the input parameters.}
#' \item{Beta_jack}{If \code{jack} is \code{TRUE}, the jackknife estimates of the regression parameters.}
#' \item{Var_jack}{If \code{jack} is \code{TRUE}, the estimated variance of the jackknife estimates.}
#' }
#' @references
#' Firth, D. and Sammut, F. (2023) Analysis of composition on the original scale of measurement. arXiv preprint arXiv:2312.10548.
#'
#' Liang, K.-Y. and Zeger, S. L. (1986) Longitudinal data analysis using generalized linear models. Biometrika 73(1), 13–22.
#'
#' Lipsitz, S. R., Laird, N. M. and Harrington, D. P. (1990) Using the jackknife to estimate the variance of regression estimators from repeated measures studies. Communications in Statistics - Theory and Methods 19(3), 821–845.
#'
#' @export
compGeeFit <- function(Y_list,
                       X_list,
                       Beta_0,
                       rho = 0,
                       Phi = NULL,
                       criterion = c("ref", "sum"),
                       jack = FALSE,
                       max.iter = 100,
                       eps = 1e-6,
                       burnin = NULL,
                       decay_factor = 1e-2,
                       verbose = 0)
{
    if (all(criterion == c("ref", "sum")) ||
        (substr("ref", 0, nchar(criterion)) == criterion))
    {
        fit.alg <- update_beta_extensive_0
    } else
    {
        if ((substr("ref", 0, nchar(criterion)) == criterion))
        {
            fit.alg <- update_beta_extensive
        } else
        {
            stop("Invalid criterion.")
        }
    }
    result <- fit.alg(
        Y_list,
        X_list,
        Beta_0,
        rho,
        Phi,
        FALSE,
        # jackknife is done later
        max.iter,
        eps,
        burnin,
        decay_factor,
        verbose
    )
    if (!jack)
        return(result)
    
    # this for the jackknife estimate and variance estimate
    NN <- length(Y_list)
    estimates <- vector("list", NN)
    for (ii in 1:NN) {
        estimates[[ii]] <- fit.alg(Y_list[-ii],
                                   X_list[-ii],
                                   result$Beta_hat,
                                   result$rho_hat,
                                   result$Phi_hat,
                                   TRUE,
                                   1) # 1-step. Following parameters are not needed
    }
    
    b_overall <- result$Beta_hat
    result$Beta_jack <- b_overall - (NN - 1) * (Reduce(`+`, estimates) / NN - b_overall)
    
    b_dot <- c(t(Reduce(`+`, estimates) / NN))
    result$Var_jack <- (NN - 1) / NN * Reduce(`+`, lapply(estimates, function(x)
        (c(t(
            x
        )) - b_dot) %*% t(c(t(
            x
        )) - b_dot)))
    return(result)
}