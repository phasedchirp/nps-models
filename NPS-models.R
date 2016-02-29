library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Load data
scores <- read.csv()

# First, get a quick estimate from the actual observations:
scores$binned <- "detractor"
scores$binned[scores$response%in%c(7,8)] <- "passive"
scores$binned[scores$response%in%c(9,10)] <- "promoter"
scores$binned <- as.factor(scores$binned)
pPro <- sum(scores$binned=="promoter")/length(scores$binned)
pDet <- sum(scores$binned=="detractor")/length(scores$binned)
pPas <- sum(scores$binned=="passive")/length(scores$binned)
(npsAll <- pPro  - pDet)


# Stan takes data as lists. This list has just responses, no grouping variable.
# The alpha value here is equivalent to adding one observation to each cell.

scoresListNoGroups <- list(N = length(scores$response),
                           K = 11,
                           y=scores$response+1, # assuming a 0-10 scale here. Drop +1 if scale starts at 1.
                           alpha = rep(2,11)) # For no smoothing, use alpha = rep(1,11)

# version with grouping variable, assuming group is a factor:
scoresListGrouped <- list(N = length(scores$response), K = 11,y=scores$response+1, alpha = rep(2,11),
                  g = as.numeric(scores$groupVar), G=length(levels(scores$groupVar)))


# Compile models:
stanModMultinomial <- stan_model(file = "NPS-multinomial.stan", model_name = "multinomial")
stanModMultinomialGrouped <- stan_model(file = "NPS-multinomial-grouped.stan", model_name = "multinomial")
stanModOrdered <- stan_model(file = "NPS-cumulative.stan", model_name = "ordered")

# Fit models:

npsMultinomial <- sampling(stanModMultinomial,data = scoresList,chains=4,iter=4000)
# Warnings on this one are *usually* harmless, but best to check convergence diagnostics here.
npsOrdered <- sampling(stanModOrdered,data = scoresList,chains=5,iter=5000)

# Compare results:
round(summary(npsMultinomial,pars="score")$summary*100,2)[,c("mean","2.5%","97.5%")]
round(summary(npsOrdered,pars="score")$summary*100,2)[,c("mean","2.5%","97.5%")]

# and compare to a convenient analytic approximation (see http://stats.stackexchange.com/questions/18603):
npsSD <- sqrt(pPro*(1-npsAll)^2 + pPas*(0-npsAll)^2 + pDet*(-1-npsAll)^2)
npsSE <- npsSD/sqrt(length(scores$binned))
round(c(100*(npsAll-npsSE),100*(npsAll+npsSE)),2)

# Looking at scores for multiple groups split on some demographic variable:
npsMultinomialGrouped <- sampling(stanModMultinomialGrouped,data = scoresList,chains=4,iter=4000)
round(summary(npsMultinomialGrouped,pars="score")$summary*100,2)[,c("mean","2.5%","97.5%")]


