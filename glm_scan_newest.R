#test by adding a space

#glm scan test with closed form solution for alternative glms

#try changing apply functions to parapply for parallel computing

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
  #why is this using ein? 
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
  
  #calculate alpha and beta for the glms
  alpha <- log((ty - yin)/off_out)
  beta <- log(yin/(ty-yin)) - log(off_in/off_out) 
  
  #calculate fitted values with alternative models
  fitted_vals <- lapply(seq(1, length(alpha)), function(i){
    exp(alpha[i] + beta[i]*indmat[,i] + log(offset))
  })
  
  #likelihood ratio calculations
  f0sum <- sum(dpois(cases, lambda = fitted(g0), log = TRUE))
  f1sums <- sapply(fitted_vals, function(object) {
    sum(dpois(cases,
              lambda = unlist(object),
              log = TRUE))
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
                          pvalue = mc, #options: mc or gumbel
                          ubpop = 0.2,
                          longlat = FALSE, cl = NULL){
  
  #ADD IN ARGUMENT CHECKING
  
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
    g0 <- glm(formula = ysim ~ PCTOWNHOME + PCTAGE65P,
              offset = log(pop), data =  data, family = family)
    
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
  
  if (pvalue == "mc"){
    pvalue <- smerc::mc.pvalue(tobs, unlist(tsim))
  }
  
  # if (pvalue == "gumbel"){
  #   #to approximate the gumbel distribution, we need the sample mean and sample
  #   #standard deviation for our generated test statistics
  #   test_stats <- unlist(tsim)
  #   
  #   #method of moments estimators
  #   mean <- mean(test_stats)
  #   sd <- sd(test_stats)
  #   beta_hat <- (sd * sqrt(6))/pi
  #   mu_hat <- mean - (0.5772 * beta_hat)
  #   
  #   #approximate distribution, right-tailed
  #   g_cdf <- function(x, mu, beta){
  #     exp(-exp(-(x-mu)/beta))
  #   }
  #   
  #   #p-value
  #   pvalue <- unlist(unname(lapply(tobs, function(x){
  #     1 - g_cdf(x = x, mu = mu_hat, beta = beta_hat)
  #   })))
  # }
  
  pruned <- smerc::sig_noc(tobs = tobs, zones = zones, pvalue = pvalue,
                           alpha = 0.1)
  
  return(pruned)
}

#test
library(smerc)
data(nypoly)

nypoly$Cases <- floor(nypoly$Cases)
#formula <- Cases ~ PCTOWNHOME + PCTAGE65P
formula <- Cases ~ 1
family <- poisson
data <- nypoly
coords <- matrix(c(nypoly$X, nypoly$Y), ncol = 2)
cases <- nypoly$Cases
pop <- nypoly$POP8
ubpop <- 0.2
pvalue <- "mc"
nsim <- 1000

test <- glm_scan_test(formula = formula, family = family, data = data, coords = coords,
                      cases = cases, pop = pop, ubpop = ubpop, pvalue = pvalue, 
                      nsim =nsim)
#only returning the MLC, not all significant ones, but it's the same as test_scan

#compare to scan.test 
test_scan <- scan.test(coords = coords, cases = cases, pop = pop, nsim = 1000, 
                       alpha = 0.1, ubpop = 0.2)


