library(deSolve)
parms   <- c(k1 = 0.04, k2 = 1e4, k3 = 3e7)
my.atol <- c(1e-6,  1e-10,  1e-6)
times <- c(0)
for(i in 1:12) {
  times[i+1] = .4 * 10 ^(i-1)
}
lsexamp <- function(t, y, p) {
  with(as.list(c(y,p)), {
    yd1 <- -k1 * y[1] + k2 * y[2]*y[3]
    yd3 <- k3 * y[2]^2
    yd2 <- -yd1-yd3
    list(c(yd1, yd2, yd3))
  })
}
out <- lsoda(c(1, 0, 0), times, lsexamp, parms, rtol = 1e-4,
             atol = my.atol)
out