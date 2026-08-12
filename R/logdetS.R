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

# .rank <- function(x) {
# ev <- eigen(x, symmetric=TRUE, only.values=TRUE)
# sum(ev$values > max(ev$values) * .Machine$double.eps^.8)
# }

# .rank_okay <- function(S, vectors = FALSE, tol = .Machine$double.eps^(2/3)) {
.rank_okay <- function(S, vectors = FALSE, tol = sqrt(.Machine$double.eps)) {
  eS <- eigen(S, symmetric = TRUE, only.values = !vectors)
  ev <- pmax(eS$values, 0)          # guard against tiny negative values
  if (ev[1] == 0) return(rep(FALSE, nrow(S)))
  out <- ev > ev[1] * tol
  if (vectors)
    attr(out, 'vectors') <- eS$vectors
  out
}

.rank <- function(S) {
  sum(.rank_okay(S))
}

.rank <- function(S) {
  qr(S)$rank
}

.logdetS_multi <- function(rho, smthi, deriv) {
  
  Si_list <- smthi$S
  
  lambda <- exp(rho)
  q <- nrow(Si_list[[1]])
  M <- length(Si_list)
  gamma <- seq_len(M)
  
  eps <- .Machine$double.eps^(2/3)

  S <- Reduce("+", Map("*", Si_list, lambda))
  
  # check rank of S
  if (.rank(S) < q) {
    Si_list_scaled <- lapply(Si_list, function(x) x / norm(x, 'F'))
    S_scaled <- Reduce("+", Si_list_scaled)
    eS_scaled <- eigen(S_scaled, symmetric = TRUE)
    # okay <- eS_scaled$values > max(eS_scaled$values) * .Machine$double.eps^.8
    # okay <- .rank_okay(S_scaled, vectors = TRUE)
    # U_positive <- attr(okay, 'vectors')[, okay, drop = FALSE]
    okay <- 1:.rank(S_scaled)#.rank_okay(S_scaled, vectors = TRUE)
    U_positive <- eS_scaled$vectors[, okay, drop = FALSE]
    Si_list <- lapply(Si_list, function(x) crossprod(U_positive, crossprod(x, U_positive)))
    q <- nrow(Si_list[[1]])
  }
  
  S <- Reduce("+", Map("*", Si_list, lambda))

  # initialize
  K <- 0
  Q <- q
  Si_tilde_list <- Si_list
  cond <- TRUE
  
  norms <- sapply(Si_tilde_list, norm, "F")
  
  while (cond) {
    # step 1
    Omega <- norms[gamma] * lambda[gamma]
    # step 2
    alpha <- gamma[Omega >= eps * max(Omega)]
    gamma_prime <- gamma[Omega < eps * max(Omega)]
    # step 3
    Si_alpha_tilde_list <- Si_tilde_list[alpha]
    Si_alpha_tilde_list_scaled <- lapply(Si_alpha_tilde_list, function(x) x / norm(x, 'F'))
    S_alpha_scaled <- Reduce("+", Si_alpha_tilde_list_scaled)
    r <- .rank(S_alpha_scaled)
    # step 4
    if (r != Q) {
      # step 5
      S_alpha <- Reduce("+", Map("*", Si_alpha_tilde_list, lambda[alpha]))
      eS_alpha <- eigen(S_alpha, symmetric = TRUE)
      U <- eS_alpha$vectors
      Ur <- U[,seq_len(r), drop = FALSE]
      Un <- U[,-seq_len(r), drop = FALSE]
      # step 6
      T_gamma_prime <- .bdiag(list(diag(K), U))
      Sp <- crossprod(T_gamma_prime, crossprod(S, T_gamma_prime))
      # step 7
      T_alpha <- .bdiag(list(diag(K), cbind(Ur, 0 * Un)))
      for (i in alpha) {
        Si_list[[i]] <- crossprod(T_alpha, crossprod(Si_list[[i]], T_alpha))
      }
      for (i in gamma_prime) {
        Si_list[[i]] <- crossprod(T_gamma_prime, crossprod(Si_list[[i]], T_gamma_prime))
        # step 8
        Si_tilde_list[[i]] <- crossprod(Un, crossprod(Si_tilde_list[[i]], Un))
      }
      # step 9
      K <- K + r
      Q <- Q - r
      S <- Sp
      gamma <- gamma_prime
    } else {
      cond <- FALSE
    }
  }
  
  qrS <- qr(S)
  out <- list(d0 = sum(log(abs(diag(qr.R(qrS))))))
  
  if (deriv > 0) {
    Si_list <- lapply(Si_list, function(x) qr.solve(qrS, x))
    out$d1 <- lambda * sapply(Si_list, function(x) sum(diag(x)))
    if (deriv > 1) {
      out$d2 <- diag(out$d1, M)
      for (i in seq_len(M)) {
        for (j in seq_len(i)) {
          out$d2[i, j] <- out$d2[i, j] - lambda[i] * lambda[j] * sum(t(Si_list[[i]]) * Si_list[[j]])
          if (i != j) 
            out$d2[j, i] <- out$d2[i, j]
        }
      }
    }
  }
  
  out
  
}



