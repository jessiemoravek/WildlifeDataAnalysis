#----------------------------------------------------
#
#  Nimble functions
#
#----------------------------------------------------

# Binomial distribution for row vectors
dbin_by_row <- nimbleFunction(
  run = function(x = double(1), P = double(1), K = double(1), log = integer(0, default = 0)) {
    J <- length(x)
    ans <- 0.0
    for(j in 1:J)
      ans <- ans + dbinom(x[j], K[j], P[j], 1)
    returnType(double())
    if(log) return(ans)
    else return(exp(ans))
  })

rbin_by_row  <- nimbleFunction(
  run = function(n = integer(), P = double(1), K = double(1)) {
    declare(ans, double(1))
    J <- length(P)
    setSize(ans, J)
    for(j in 1:J)
      ans[j] <- rbinom(1, K[j], P[j])
    returnType(double(1))
    return(ans)
  })

registerDistributions(list(
  dbin_by_row = list(
    BUGSdist = "dbin_by_row(P, K)",
    Rdist = "dbin_by_row(P, K)",
    range = c(0, Inf),
    types = c('value = double(1)', 'P = double(1)', 'K = double(1)'))
))
#-------------------------------
# Poisson distribution

dpois_by_row <- nimbleFunction(
  run = function(x = double(1), bigLambda = double(1), K = double(1), log = integer(0, default = 0)) {
    J <- length(x)
    ans <- 0.0
    for(j in 1:J)
      ans <- ans + dpois(x[j], bigLambda[j]*K[j], 1)
    returnType(double())
    if(log) return(ans)
    else return(exp(ans))
  })

rpois_by_row  <- nimbleFunction(
  run = function(n = integer(), bigLambda = double(1), K = double(1)) {
    declare(ans, double(1))
    J<- length(bigLambda)
    setSize(ans, J)
    for(j in 1:J)
      ans[j] <- rpois(1, bigLambda[j]*K[j])
    returnType(double(1))
    return(ans)
  })

registerDistributions(list(
  dpois_by_row = list(
    BUGSdist = "dpois_by_row(bigLambda, K)",
    Rdist = "dpois_by_row(bigLambda, K)",
    range = c(0, Inf),
    types = c('value = double(1)', 'bigLambda = double(1)', 'K = double(1)'))
))

#---------------------------------------------------
# for discrete state space. input is an sf object
make_grid <- function(x, cell_diameter, cell_area, square= FALSE, clip = FALSE, offset=c(0,0)) {
  if (missing(cell_diameter)) {
    if (missing(cell_area)) {
      stop("Must provide cell_diameter or cell_area")
    } else {
      cell_diameter <- sqrt(2 * cell_area / sqrt(3))
    }
  }
  # generate array of hexagon centers
  g <- st_make_grid(x, cellsize = cell_diameter, what="polygons", square=square, offset = offset)
  
  # clip to boundary of study area
  if (clip) {
    g <- st_intersection(g, x)
  }
  g<- st_sf(HexID = 1:length(g), geometry=g)
  return(g)
}
