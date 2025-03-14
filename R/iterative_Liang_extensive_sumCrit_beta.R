update_beta_extensive <- function(Y_list,
                                  X_list,
                                  Beta_new,
                                  rho,
                                  Phi = NULL,
                                  jack = FALSE,
                                  max.iter = 100,
                                  eps = 1e-6,
                                  burnin = NULL,
                                  lambda_0 = 1e-2,
                                  verbose = 0)
{
    # Y_list is the list of time series, each of which is a Ti x D matrix
    # X_list is the list of time series covariates, each of which is a Ti x (P+1) matrix -- MM = P dimensions
    NN <- length(Y_list)
    TT <- sapply(Y_list, function(yy)
        dim(yy)[1])
    DD <- dim(Y_list[[1]])[2]
    MM <- dim(X_list[[1]])[2] # this is P+1, I call it MM because the variable PP has other uses in this code
    CC <- diag(1, DD) - 1 / DD
    if (is.null(MM))
        MM <- 1
    
    # other initializations
    if (is.null(Phi))
    {
        Phi <- diag(1, DD)
    } else{
        if (length(Phi) == 1)
        {
            Phi <- diag(Phi, DD)
        }
    }
    if (is.null(burnin))
        burnin = max.iter
    
    # initialize stop conditions
    iter <- 0
    diff <- 0
    lambda <- NA
    out_diff <- NULL
    conv <- FALSE
    
    # fitted and residuals for correlations
    fit_list <- vector("list", NN)
    res_list <- vector("list", NN)

    starttime <- Sys.time()
    while (iter < max.iter & !conv) {
        if (verbose && (iter %% verbose) == 0)
        {
            cat("\n\nIteration: ", iter)
            cat("\nDifference: ", diff)
            cat("\nLambda: ", lambda)
            cat("\nBeta:\n")
            print(Beta_new)
            cat("\nRho:", rho)
            cat("\nPhi:\n")
            print(Phi)
        }
        iter <- iter + 1
        
        # initialize sums
        sd <- 0
        ss <- 0
        svar <- 0
        Beta_curr <- Beta_new
        
        # cicle for update beta
        for (ii in 1:NN) {
            # initialize matrices/vectors
            Ti <- TT[ii]
            X_i <- X_list[[ii]]
            tot_yi <- rowSums(Y_list[[ii]]) # totals for residuals ecc
            yi <- c(t(Y_list[[ii]])) # vec for beta update
            XXi <- NULL
            vec_pi <- NULL
            PP <- NULL
            VV <- NULL
            res_mat <- NULL
            
            for (jj in 1:Ti) {
                # covariates at time jj
                xx <- X_i[jj, ]
                # large user covariate matrix
                XXi <- rbind(XXi, kronecker(diag(1, DD), t(xx)))
                
                # eta, pi at time jj
                eta <- apply(Beta_curr, 1, function(bb)
                    sum(xx * bb))
                pi <- softmax(eta)
                # user expected vector / all_users fit matrix for correlation
                vec_pi <- c(vec_pi, tot_yi[jj] * pi) # extensive modification here
                
                # multinomial covariance matrix
                Pi <- diag(pi) - pi %*% t(pi)
                if (jj == 1)
                {
                    # initialization
                    PP <- Pi
                    VV <- Pi %*% Phi %*% Pi
                } else
                {
                    # large user multinomial covariance matrix
                    PP <- Matrix::bdiag(PP, Pi)
                    # generalized Wedderburn variance function
                    VV <- Matrix::bdiag(VV, Pi %*% Phi %*% Pi)
                }
            }
            # build V tilde
            RR <- kronecker(rho^abs(matrix(1:Ti, Ti, Ti, TRUE) - 1:Ti), diag(1, DD))
            eig_VV <- eigen(VV)
            sq_VV <- eig_VV$vectors %*% sqrt(diag(abs(eig_VV$values))) %*% solve(eig_VV$vectors) # square of V
            V_tilde <- Re(sq_VV %*% RR %*% sq_VV)
            
            # get Sd, S
            dtvm <- t(XXi) %*% PP %*% MASS::ginv(V_tilde)
            sd <- sd + dtvm %*% PP %*% XXi
            ss <- ss + dtvm %*% (yi - vec_pi)
        }
        Beta_curr <- c(t(Beta_curr))
        upd <- MASS::ginv(as.matrix(sd)) %*% ss
        
        lambda <- max(1 / max(Matrix::norm(upd, "M"), 1), 1e-6)
        # decay
        if (iter > burnin)
            lambda <- lambda * exp(-lambda_0 * (iter - burnin))
        Beta_new <- Beta_curr + lambda * upd
        
        # diff <- norm(upd, "e") / norm(Beta_new, "e")
        diff <- max(Beta_new - Beta_curr)
        out_diff <- c(out_diff, diff)
        
        Beta_new <- matrix(Beta_new, DD, byrow = TRUE)
        if(jack)
            return(Beta_new)
        # cycle for update phi, rho
        for (ii in 1:NN) {
            # initialize matrices/vectors
            Ti <- TT[ii]
            X_i <- X_list[[ii]]
            tot_yi <- rowSums(Y_list[[ii]])
            yi <- c(t(Y_list[[ii]] / tot_yi)) # de-extensive
            mat_pi <- NULL
            res_mat <- NULL
            for (jj in 1:Ti) {
                # covariates at time jj
                xx <- X_i[jj, ]
                
                # eta, pi at time jj
                eta <- apply(Beta_new, 1, function(bb)
                    sum(xx * bb))
                pi <- softmax(eta)
                # user expected vector / all_users fit matrix for correlation
                mat_pi <- rbind(mat_pi, pi) # de-extensive modification
                # multinomial covariance matrix
                Pi <- diag(pi) - pi %*% t(pi)
                
                res_mat <- rbind(res_mat, t(MASS::ginv(Pi) %*% (yi[((jj - 1) * DD + 1):(jj * DD)] - pi))) # de-extensive modification
                
            }
            fit_list[[ii]] <- mat_pi
            res_list[[ii]] <- res_mat
        }
        
        # correlation based on residuals
        past_fit <- c(unlist(sapply(res_list, function(x)
            c(x[-dim(x)[1], ]))))
        future_fit <- c(unlist(sapply(res_list, function(x)
            c(x[-1, ]))))
        rho <- stats::cor(past_fit, future_fit)
        
        # update phi
        Phi <- Reduce(`+`, lapply(res_list, function(res_mat)
            Reduce(
                `+`,
                apply(res_mat, 1, function(res_row)
                    res_row %*% t(res_row), simplify = FALSE)
            ))) / (sum(TT) - (DD - 1) * MM) # gradi di libertà: P * DD - 1
        
        if (diff < eps)
        {
            conv <- TRUE
        }
    }
    
    alg.params <- list(max.iter = max.iter,
                       burnin = burnin,
                       eps = eps)
    return(
        list(
            Beta_hat = Beta_new,
            rho_hat = rho,
            Phi_hat = Phi,
            conv = conv,
            diff = diff,
            out_diff = out_diff,
            n_iter = iter,
            time = Sys.time() - starttime,
            alg.params = alg.params
        )
    )
}
