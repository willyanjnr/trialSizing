## ============================================================================
## data-raw/uniformity_trial.R
## Generates the simulated example dataset `uniformity_trial`.
##
## A uniformity trial is a field sown uniformly (one genotype, one management)
## and harvested in a fine grid of small basic experimental units (BEU). The
## spatial variation left over is pure "environmental" noise, and it is what
## plot-size methods use to decide how large a plot must be.
##
## This dataset is fully SIMULATED for teaching, so it can be redistributed
## freely. Each of the three trials is an 8 x 12 grid of 1 m^2 BEU whose yield
## (g m^-2) is drawn from a separable first-order autoregressive (AR(1)) field
## -- a spatially correlated "patch" component that persists when neighbouring
## units are pooled (this produces the CV plateau) plus an independent nugget
## that averages out quickly (this produces the steep early decline). The
## parameters were tuned so the CV(x) curve looks like a real chickpea/oat
## uniformity trial: CV ~ 23% at x = 1 m^2, a plateau near 6-8%, and the usual
## MCM < LRP < QRP ordering of the optimum.
##
## Run with:  Rscript --vanilla data-raw/uniformity_trial.R
## ============================================================================

sim_grid <- function(nr, nc, mu, rho, sd_patch, sd_nugget) {
  ar1 <- function(n, rho) {
    z <- numeric(n); z[1] <- rnorm(1)
    for (i in 2:n) z[i] <- rho * z[i - 1] + sqrt(1 - rho^2) * rnorm(1)
    z
  }
  R <- ar1(nr, rho); C <- ar1(nc, rho)
  patch <- outer(R, C, "+") / sqrt(2)
  patch <- (patch - mean(patch)) / stats::sd(as.vector(patch))
  nug <- matrix(rnorm(nr * nc), nr, nc)
  round(mu + sd_patch * patch + sd_nugget * nug, 2)
}

set.seed(101)
grids <- lapply(1:3, function(i)
  sim_grid(nr = 8, nc = 12, mu = 250, rho = 0.75,
           sd_patch = 26, sd_nugget = 58))
names(grids) <- paste0("T", 1:3)

## Wide "field map" layout: one row per grid row, one column per grid column.
uniformity_trial <- do.call(rbind, lapply(names(grids), function(nm) {
  m <- grids[[nm]]
  cols <- as.data.frame(m)
  names(cols) <- sprintf("col%02d", seq_len(ncol(m)))
  cbind(data.frame(trial = nm, row = seq_len(nrow(m)),
                   stringsAsFactors = FALSE), cols)
}))
rownames(uniformity_trial) <- NULL

save(uniformity_trial, file = "data/uniformity_trial.rda", compress = "xz")
