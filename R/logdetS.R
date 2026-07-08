## functions for calculating log(|S|)

.logdetS <- function(smth, rho, deriv=0) {
id <- rep(seq_along(smth), sapply(smth, function(x) length(x$S)))
rho <- split(rho, id)
out <- ldSi <- list()
for (i in seq_along(rho)) {
  if (length(rho[[i]]) == 1) {
    ldSi[[i]] <- .logdetS_single(rho[[i]], smth[[i]], deriv)
  } else 
    ldSi[[i]] <- .logdetS_multi(rho[[i]], smth[[i]], deriv)
}
out$d0 <- Reduce("+", lapply(ldSi, function(x) x$d0))
if (deriv > 0) {
  out$d1 <- unlist(lapply(ldSi, function(x) x$d1))
  if (deriv > 1) 
    out$d2 <- .bdiag(lapply(ldSi, function(x) x$d2))
}
out
}

.logdetS_single <- function(rho, smthi, deriv) {
out <- list(d0 = rho * smthi$rank)
if (deriv > 0) {
  out$d1 <- smthi$rank
  if (deriv > 1) 
    out$d2 <- matrix(0, length(rho), length(rho))
}
out
}

.rank <- function(x) {
ev <- eigen(x, symmetric=TRUE, only.values=TRUE)
sum(ev$values > max(ev$values) * .Machine$double.eps^.8)
}

# .logdetS_multi <- function(rho, smthi, deriv) {
# 
# Sl <- smthi$S
# 
# # this calculate log(|S|_+) using similarity transforms as described in
# # Appendix B of Wood, S. N. (JRSSB, 2011). "Fast stable REML and MLE of semiparametric GLMs"
# 
# lambda <- exp(rho)
# q <- nrow(Sl[[1]])
# eps <- .Machine$double.eps^(1/3)
# 
# S <- Reduce("+", Map("*", Sl, lambda))
# 
# # check rank of S
# if (.rank(S) < q) {
#   tS <- Sl
#   for (i in seq_along(tS))
#     tS[[i]] <- tS[[i]] / norm(tS[[i]], "F")
#   S <- Reduce("+", Map("*", tS, lambda))
#   eS <- eigen(S, symmetric=TRUE)
#   okay <- eS$values > max(eS$values) * .Machine$double.eps^.8
#   Up <- eS$vectors[, okay, drop=FALSE]
#   Sl <- lapply(Sl, function(x) crossprod(Up, crossprod(x, Up)))
#   q <- nrow(Sl[[1]])
# }
# 
# S <- Reduce("+", Map("*", Sl, lambda))
# M <- length(Sl)
# 
# # initialize
# Q <- q
# K <- 0
# gamma <- seq_len(M)
# Sb <- Sl
# cond <- TRUE
# 
# while (cond) {
#   # step 1
#   norms <- sapply(Sb[gamma], norm, "F")
#   Omega <- try(norms * lambda[gamma])
#   if (inherits(Omega, "try-error")) browser()
#   # step 2
#   alpha <- gamma[Omega >= eps * max(Omega)]
#   gammap <- gamma[Omega < eps * max(Omega)]
#   # step 3
#   Sa <- Reduce("+", Map("*", Sb[alpha], lambda[alpha]))
#   r <- .rank(Sa)
#   # step 4
#   if (r != Q) {
#     # step 5
#     eS <- eigen(Sa, symmetric=TRUE)
#     U <- eS$vectors
#     Ur <- U[,seq_len(r), drop=FALSE]
#     Un <- U[,-seq_len(r), drop=FALSE]
#     # step 6
#     Tg <- .bdiag(list(diag(K), U))
#     Sp <- try(crossprod(Tg, crossprod(S, Tg)), silent= TRUE)
#     # step 7
#     Ta <- cbind(.bdiag(list(diag(K), Ur)), 0 * Un)
#     for (i in alpha)
#       Sl[[i]] <- crossprod(Ta, crossprod(Sl[[i]], Ta))
#     for (i in gammap)
#       Sl[[i]] <- crossprod(Tg, crossprod(Sl[[i]], Tg))
#     # step 8
#     for (i in gammap)
#       Sb[[i]] <- crossprod(Un, crossprod(Sb[[i]], Un))
#     # step 9
#     K <- K + r
#     Q <- Q - r
#     S <- Sp
#     gamma <- gammap
#   } else {
#     cond <- FALSE
#   }
# }
# 
# cS <- chol(S)
# out <- list(d0 = 2 * sum(log(diag(cS))))
# 
# if (deriv > 0) {
#   iS <- chol2inv(cS)
#   if (deriv > 1) {
#     Sl <- lapply(Sl, function(x) crossprod(iS, x))
#     out$d1 <- lambda * sapply(Sl, function(x) sum(diag(x)))
#     out$d2 <- diag(out$d1, M)
#     for (i in seq_len(M)) {
#       for (j in seq_len(i)) {
#         out$d2[i, j] <- out$d2[i, j] - lambda[i] * lambda[j] * sum(t(Sl[[i]]) * Sl[[j]])
#         if (i != j) 
#           out$d2[j, i] <- out$d2[i, j]
#       }
#     }
#   } else {
#     out$d1 <- lambda * sapply(Sl, function(x) sum(iS * x))
#   }
# }
# 
# out
# 
# }

.logdetS_multi <- function(rho, smthi, deriv = 0, tol = .Machine$double.eps^(1/3)) {
  # rho   : vector of length m containing the log-smoothing parameters (rho_j)
  # Sl    : a list of m symmetric positive semi-definite penalty matrices (S_j)
  # deriv : 0 for log-det only, 1 for gradient, 2 for gradient + Hessian
  
  Sl <- smthi$S
  
  if (!deriv %in% c(0, 1, 2)) {
    stop("Argument 'deriv' must be 0, 1, or 2.")
  }
  
  lambda <- exp(rho)
  m <- length(lambda)
  p <- nrow(Sl[[1]])
  
  # Step 1: Pre-transform each individual block to isolate its unique range space
  active_dims <- integer(0)
  Sl_transformed <- list()
  
  for (j in 1:m) {
    S_j <- Sl[[j]]
    ev <- eigen(S_j, symmetric = TRUE)
    rank_j <- sum(ev$values > (max(ev$values) * tol))
    
    if (rank_j > 0) {
      U_j <- ev$vectors
      S_j_prime <- t(U_j) %*% S_j %*% U_j
      Sl_transformed[[j]] <- S_j_prime
      active_dims <- union(active_dims, 1:rank_j)
    } else {
      Sl_transformed[[j]] <- matrix(0, p, p)
    }
  }
  
  # Initialize results structures
  d0 <- -Inf
  d1 <- if (deriv >= 1) numeric(m) else NULL
  d2 <- if (deriv == 2) matrix(0, m, m) else NULL
  
  # Step 2: Combine the pre-stabilised blocks
  if (length(active_dims) > 0) {
    S_lambda_active <- matrix(0, length(active_dims), length(active_dims))
    for (j in 1:m) {
      S_lambda_active <- S_lambda_active + lambda[j] * Sl_transformed[[j]][active_dims, active_dims, drop = FALSE]
    }
    
    # Step 3: Conditional Eigen-decomposition to save calculations
    need_vectors <- (deriv >= 1)
    ev_final <- eigen(S_lambda_active, symmetric = TRUE, only.values = !need_vectors)
    
    pos_idx <- which(ev_final$values > (max(ev_final$values) * tol))
    
    if (length(pos_idx) > 0) {
      positive_evs <- ev_final$values[pos_idx]
      d0 <- sum(log(positive_evs))
      
      # Stop early if no derivatives are requested
      if (deriv == 0) {
        return(list(d0 = d0, d1 = d1, d2 = d2))
      }
      
      # Core structures needed for derivatives
      V_pos <- ev_final$vectors[, pos_idx, drop = FALSE]
      inv_evs <- 1 / positive_evs
      
      # --- DERIVATIVE LEVEL 1 (Gradient Only) ---
      if (deriv == 1) {
        for (k in 1:m) {
          S_k_active <- Sl_transformed[[k]][active_dims, active_dims, drop = FALSE]
          # Vectorized diagonal tracking avoiding full matrix products
          diag_M <- colSums(V_pos * (S_k_active %*% V_pos))
          d1[k] <- lambda[k] * sum(inv_evs * diag_M)
        }
      }
      
      # --- DERIVATIVE LEVEL 2 (Gradient + Hessian Matrix) ---
      if (deriv == 2) {
        # Precompute the fully projected M_k matrices to avoid repeating O(n^3) ops
        M_list <- list()
        for (k in 1:m) {
          S_k_active <- Sl_transformed[[k]][active_dims, active_dims, drop = FALSE]
          M_list[[k]] <- t(V_pos) %*% S_k_active %*% V_pos
          
          # Calculate the gradient term directly from the diagonal of M_k
          d1[k] <- lambda[k] * sum(inv_evs * diag(M_list[[k]]))
        }
        
        # Populate the symmetric Hessian matrix
        # d^2 / (d rho_l d rho_k) = delta_kl * d1_k - lambda_k * lambda_l * tr(S^-1 * S_l * S^-1 * S_k)
        for (k in 1:m) {
          for (l in k:m) {
            # Matrix trace expression expressed efficiently via element-wise matrix multiplication
            # tr(diag(inv_evs) %*% M_l %*% diag(inv_evs) %*% M_k)
            # which simplifies to: sum_{i,j} (inv_evs[i] * inv_evs[j] * M_l[i,j] * M_k[i,j])
            trace_term <- sum((inv_evs %*% t(inv_evs)) * M_list[[l]] * M_list[[k]])
            
            val <- -lambda[k] * lambda[l] * trace_term
            if (k == l) {
              val <- val + d1[k]
            }
            
            d2[k, l] <- val
            if (k != l) {
              d2[l, k] <- val  # Enforce symmetry mapping
            }
          }
        }
      }
    }
  }
  
  return(list(
    d0 = d0,  # Log-determinant
    d1 = d1,  # Gradient vector (w.r.t rho) or NULL
    d2 = d2   # Hessian matrix (w.r.t rho) or NULL
  ))
}

