#quasi poisson glm function


#function for creating the matrix indicating which regions are in which zones
#number of columns = number of zones, number of rows = number of regions
zones2ind <- function(zones, nregions) {
  ilist <- lapply(zones, function(x) {
    u <- numeric(nregions)
    u[x] <- 1
    u
  })
  imat <- do.call(cbind, ilist)
  unname(imat)
}


#function for finding the MLCs
glm_scan_stat <- function(cases,
                          ty,
                          yin,
                          ein,
                          pop,
                          nn,
                          indmat,
                          g0) {
  
  #setup
  
  #indicator for test statistic
  tstat_ind <- yin/ein > (sum(cases) - yin)/(sum(cases) - ein)
  
  #alternative glms use the same beta*X (beta from g0) (as an additional offset)
  #and include an indicator for candidate zone
  
  #if there are no covariates, there is no additional offset
  if (dim(g0$model)[2] == 2){
    offset <- pop
  } else { 
    #calculate what the alternative offset is when there are covariates
    X <- data.matrix(g0$model[2:(length(g0$model)-1)])
    coef <- matrix(g0$coefficients[-1], ncol = 1)
    alt_offset <- X %*% coef
    offset <- exp(alt_offset + log(pop))
  }
  
  #calculate the total offset in and outside of each zone
  off_in <- smerc::nn.cumsum(nn, offset)
  off_out <- sum(offset) - off_in
  
  #calculate alpha and beta for the alt glms
  alpha <- log((ty - yin)/off_out)
  beta <- log(yin/(ty-yin)) - log(off_in/off_out) 
  
  #calculate fitted values for alternative models
  fitted_vals <- lapply(seq(1, length(alpha)), function(i){
    exp(alpha[i] + beta[i]*indmat[,i] + log(offset))
  })
  
  
  #likelihood ratio calculations
  
  #quasi-poisson likelihood function, this is undefined if there aren't any cases
  #in a region because of the log(y)
  qp <- function(y, mu, phi){
    if (min(mu) < 0) { 
      stop("mu be have non-negative values")
    }
    tstat <- numeric(length(y)) #create vector of 0s
    yind <- (y > 0) # tstat is 0 if y == 0, so only compute for y > 0
    tstat[yind] <- 1/phi * (y[yind] * log(mu[yind]) - mu - y[yind] * log(y[yind]) + y[yind])
    # 1/phi * (y*log(mu) - mu - y*log(y) + y)
    return(tstat)
  }
  
  #indicator for which regions have cases
  ind <- (cases > 0)
  
  #likelihood under the null hypothesis
  phi0 <- sigma(g0)**2 #dispersion parameter for the null model
  
  f0sum <- sum(unname(qp(y = cases, mu = fitted(g0), phi = phi0)[ind]))
  
  #likelihood under the alternative hypotheses
  f1sums <- sapply(fitted_vals, function(object) {
    #dispersion parameter for the model
    phi_alt <- (1/(length(cases)-2))*sum(((cases - unlist(object))**2)/unlist(object))

    sum(unname(qp(y = cases, mu = unlist(object), phi = phi_alt)[ind]))
  })
  
  tobs = (f1sums - f0sum) * tstat_ind 
  
  return(tobs)
}


glm_scan_test <- function(formula,
                          family = poisson,
                          data,
                          coords,
                          cases,
                          pop,
                          nsim = 499,
                          alpha = 0.1,
                          ubpop = 0.2,
                          longlat = FALSE, cl = NULL){
  
  
  #setup
  cases <- floor(cases)
  
  #convert coords to matrix
  coords <- as.matrix(coords)
  
  #calculate distances
  d <- as.matrix(dist(coords))
  
  #create nearest neighbor list
  nn <- smerc::nnpop(d, pop = pop, ubpop = ubpop)
  
  #zones listed out, instead of being given as nested lists for each candidate zone
  zones <- smerc::nn2zones(nn)
  
  #use indmat function to create a matrix of indicator variables for each 
  #candidate zone, each column represents a zone
  indmat <- zones2ind(zones, nregions = length(nn))
  
  #calculate null glm here and feed it to the tobs function instead of calculating
  #it in both functions
  g0 <- glm(formula = formula,
            offset = log(pop), data =  data, family = family)
  
  #calculate things for glm_scan_stat function
  ex <- unname(g0$fitted.values)
  ein <- smerc::nn.cumsum(nn, ex)
  ty <- sum(cases) #sum of cases in the study area
  yin <- smerc::nn.cumsum(nn, cases)
  
  tobs <- glm_scan_stat(cases = cases, ty = ty, yin = yin,
                        ein = ein, pop = pop, nn = nn, indmat = indmat,
                        g0 = g0)
  
  #p-value calculations
  #tsim calculation
  
  tsim <- pbapply::pblapply(seq_len(nsim), function(i){
    
    #generate data set
    ysim <- stats::rmultinom(1, size = ty, prob = ex)
    yin <- smerc::nn.cumsum(nn, ysim)
    
    #calculate g0 for the simulated data set
    g0 <- glm(formula = ysim ~ pctownhome + pctage65p,
              offset = log(pop8), data =  data, family = family)
    
    ex <- unname(g0$fitted.values)
    ein <- smerc::nn.cumsum(nn, ex)
    
    t <- glm_scan_stat(cases = ysim, ty = sum(ysim), yin = yin,
                       ein = ein, pop = pop, nn = nn, indmat = indmat, 
                       g0 = g0)
    #return the max test stat
    max(t, na.rm = TRUE) #add in na.rm = TRUE, because there are zones with 
    #zero cases and this returns NA for the test stat
  })
  
  
  #detecting multiple clusters
  wdup <- nndup(nn)
  w0 <- which(tobs == 0 | yin == 0 | wdup)
  zones <- nn2zones(nn)
  zones <- zones[-w0]
  tobs <- tobs[-w0]
  
  pvalue <- smerc::mc.pvalue(tobs, unlist(tsim))
  
  pruned <- smerc::sig_noc(tobs = tobs, zones = zones, pvalue = pvalue,
                           alpha = 0.1)
  
  return(pruned)
}

#test

library(smerc)
data(nysf)

nysf$cases <- floor(nysf$cases)
formula <- cases ~ pctownhome + pctage65p
family <- quasipoisson
data <- nysf
coords <- matrix(c(nysf$x, nysf$y), ncol = 2)
cases <- nysf$cases
pop <- nysf$pop8
ubpop <- 0.2
nsim <- 100

test <- glm_scan_test(formula = formula, family = family, data = data, coords = coords,
                      cases = cases, pop = pop, ubpop = ubpop, nsim =nsim)





