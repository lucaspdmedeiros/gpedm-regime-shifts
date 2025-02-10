# Function that uses QR decomposition to compute the Lyapunov Exponents
# of a discrete-time dynamical system from a sequence of Jacobian matrices

# Arguments:
# J: a list of sequential Jacobian matrices 

# Output:
# LEs: Lyapunov Exponents

le_qr_decomp <- function(J) {
  n_sp <- nrow(J[[1]])
  n <- length(J)
  Q <- diag(rep(1, n_sp))
  R_diag <- matrix(data = NA, nrow = n, ncol = n_sp)
  for (i in 1:n){
    Phi <- J[[i]]
    QR_decomp <- qr(Phi %*% Q)
    Q <- qr.Q(QR_decomp)
    R <- qr.R(QR_decomp)
    sgn <- sign(diag(R))
    Q <- Q %*% diag(sgn)
    R <- diag(sgn) %*% R
    R_diag[i, ] <- log(diag(R))
  }
  LEs <- apply(R_diag, 2, sum) / n
  return(LEs)
}
