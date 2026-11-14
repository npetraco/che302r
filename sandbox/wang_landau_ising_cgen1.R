#Claude:
# Wang-Landau sampling for the 2D Ising model (J = 1, kB = 1, periodic BCs)
wang_landau_ising <- function(L = 8, ln_f_final = 1e-5, flat = 0.8,
                              check_every = 10000, seed = 1) {
  set.seed(seed)
  N <- L * L
  s <- matrix(sample(c(-1, 1), N, replace = TRUE), L, L)

  # Total energy: count each bond once (right and down neighbours)
  energy <- function(s) -sum(s * (s[, c(2:L, 1)] + s[c(2:L, 1), ]))
  E <- energy(s)

  # Energies are -2N, -2N+4, ..., 2N  ->  N+1 bins
  idx <- function(E) (E + 2 * N) / 4 + 1
  nb  <- N + 1
  lng <- numeric(nb)           # ln g(E), initialised to 0 (g = 1)
  H   <- numeric(nb)           # histogram
  visited <- logical(nb)       # bins actually reachable (E = -2N+4 and 2N-4 are not)
  ln_f <- 1                    # ln of the modification factor

  prev <- function(k) if (k == 1) L else k - 1
  nxt  <- function(k) if (k == L) 1 else k + 1

  while (ln_f > ln_f_final) {
    repeat {
      for (k in seq_len(check_every)) {
        # Propose a single spin flip
        i <- sample.int(L, 1); j <- sample.int(L, 1)
        dE <- 2 * s[i, j] * (s[prev(i), j] + s[nxt(i), j] +
                               s[i, prev(j)] + s[i, nxt(j)])
        a <- idx(E); b <- idx(E + dE)

        # Accept with probability min(1, g(E1)/g(E2))
        if (log(runif(1)) < lng[a] - lng[b]) {
          s[i, j] <- -s[i, j]
          E <- E + dE
        }

        # Update the CURRENT state's bin, whether or not the move was accepted
        a <- idx(E)
        lng[a] <- lng[a] + ln_f
        H[a]   <- H[a] + 1
        visited[a] <- TRUE
      }
      # Flatness: every reachable bin at least `flat` * mean(H)
      if (all(H[visited] > flat * mean(H[visited]))) break
    }
    cat(sprintf("ln f = %.2e  flat\n", ln_f))
    H[] <- 0
    ln_f <- ln_f / 2
  }

  # Normalise so that sum_E g(E) = 2^N
  lng[!visited] <- -Inf
  m <- max(lng[visited])
  lse <- m + log(sum(exp(lng[visited] - m)))
  lng <- lng + N * log(2) - lse

  list(E = -2 * N + 4 * (0:N), lng = lng, N = N)
}


# Thermodynamics from ln g(E) at any temperature
thermo <- function(res, T) {
  ok <- is.finite(res$lng)
  E <- res$E[ok]; lng <- res$lng[ok]
  sapply(T, function(t) {
    w <- lng - E / t
    w <- exp(w - max(w)); w <- w / sum(w)   # Boltzmann weights
    U  <- sum(w * E)
    C  <- (sum(w * E^2) - U^2) / t^2
    c(T = t, U_per_spin = U / res$N, C_per_spin = C / res$N)
  })
}

# ---- Run ----
res <- wang_landau_ising(L = 8)

cat("\nSanity check: g(ground state) =", exp(res$lng[1]), "(exact: 2)\n")

T <- seq(1, 5, by = 0.05)
out <- as.data.frame(t(thermo(res, T)))

# Exact Tc for the infinite lattice is 2/ln(1+sqrt(2)) ~ 2.269
op <- par(mfrow = c(1, 3))
plot(res$E, res$lng, type = "b", xlab = "E", ylab = "ln g(E)")
plot(out$T, out$U_per_spin, type = "l", xlab = "T", ylab = "U / N")
plot(out$T, out$C_per_spin, type = "l", xlab = "T", ylab = "C / N")
abline(v = 2 / log(1 + sqrt(2)), lty = 2)
par(op)
