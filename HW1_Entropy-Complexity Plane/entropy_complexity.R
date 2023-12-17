shenons_entropy <- function(P) {
    P <- P[P != 0]
    return(-sum(P * log(P)))
}

entropic_measure <- function(P) {
    S_max <- shenons_entropy(rep(1, length(P)) / length(P))
    return(shenons_entropy(P) / S_max)
}

jensen_shannon_divergence <- function(P_1, P_2) {
    N <- length(P_1)
    Q_0 <- -2 / (((N + 1) / N * log(N + 1)) - 2 * log(2 * N) + log(N))
    return(Q_0 * (shenons_entropy((P_1 + P_2) / 2) - shenons_entropy(P_1) / 2 - shenons_entropy(P_2) / 2))
}

complexity <- function(P) {
    return(jensen_shannon_divergence(P, rep(1, length(P)) / length(P)) * entropic_measure(P))
}

min_borders <- function(d, n_steps = 500) {
    N <- factorial(d)
    d_step <- (1 - 1 / N) / n_steps

    p_min <- seq(1 / N, 1, d_step)
    min_complexity <- numeric(length(p_min))
    min_entropy <- numeric(length(p_min))

    for (i in seq_along(p_min)) {
        P <- c(p_min[i], rep((1 - p_min[i]) / (N - 1), N - 1))
        min_entropy[i] <- entropic_measure(P)
        min_complexity[i] <- complexity(P)
    }

    return(list(min_entropy, min_complexity))
}

max_borders <- function(d, n_steps = 500) {
    N <- factorial(d)
    d_step <- (1 - 1 / N) / n_steps

    max_complexity <- numeric(0)
    max_entropy <- numeric(0)

    for (i in seq_len(N - 1)) {
        p_max <- seq(0, 1 / (N - i), d_step)

        for (j in seq_along(p_max)) {
            P <- c(p_max[j], rep((1 - p_max[j]) / (N - i - 1), N - i - 1))

            if (length(P) != N) {
                P <- c(P, rep(0, i))
            }

            max_entropy <- c(max_entropy, entropic_measure(P))
            max_complexity <- c(max_complexity, complexity(P))
        }
    }

    return(list(max_entropy, max_complexity))
}

ordinal_pattern <- function(x, dim) {
    ordinal_numbers <- seq(0, (dim - 1), by = 1)
    possible_pattern <- (combinat::permn(ordinal_numbers))
    result <- 0
    result[1:length(possible_pattern)] <- 0

    for (i in 1:(length(x) - (dim - 1))) {
        temp <- x[i:(i + (dim - 1))]
        tempseq <- seq(0, dim - 1, by = 1)
        tempdata <- data.frame(temp, tempseq)
        tempdata <- tempdata[order(temp), ]

        for (j in 1:length(possible_pattern)) {
            if (all(possible_pattern[[j]] == tempdata$tempseq)) {
                result[j] <- result[j] + 1
            }
        }
    }

    return(result)
}


# Skew tent map
skew_tent_map <- function(n = 2^15, omega = 0.1847, x0 = 0.5) {
  x <- numeric(n)
  x[1] <- x0
  for (i in 2:n) {
    if (x[i-1] < omega) {
      x[i] <- x[i-1] / omega
    } else {
      x[i] <- (1 - x[i-1]) / (1 - omega)
    }
  }
  return(x)
}

# Logistic map
logistic_map <- function(n = 2^15, r = 4, x0 = 0.4) {
  x <- numeric(n)
  x[1] <- x0
  for (i in 2:n) {
    x[i] <- r * x[i-1] * (1 - x[i-1])
  }
  return(x)
}

# Schuster's map
schuster_map <- function(n = 2^15, z = 2, x0 = 0.5) {
  x <- numeric(n)
  x[1] <- x0
  for (i in 2:n) {
    x[i] <- x[i-1] + (x[i-1]^z)
    x[i] <- x[i] %% 1
  }
  return(x)
}

# Henon map
henon_map <- function(n = 1000000, a = 1.4, b = 0.3, x0 = 0.4) {
  x <- numeric(n)
  x[1] <- x0
  for (i in 2:n) {
    x[i] <- 1 - a * x[i-1]^2 + b * x[i-1]
  }
  return(x)
}

sine_data <- sin(seq(0, 1000, 0.01))
gaussian_noise <- rnorm(100000)


min_data <- min_borders(6)
X_min <- min_data[[1]]
y_min <- min_data[[2]]

max_data <- max_borders(6)
X_max <- max_data[[1]]
y_max <- max_data[[2]]

plot(X_max, y_max, type = "l", col = "blue", xlab = "X", ylab = "Y", main = "Entropy/complexity plane")
lines(X_min, y_min, col = "blue")
p_sine <- ordinal_pattern(sine_data)
points(entropic_measure(p_sine), complexity(p_sine), col = "red")
print(entropic_measure(p_sine))
