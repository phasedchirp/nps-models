data {
    int<lower=1> N; // number of observations
    int<lower=2> K; // number of bins
    vector[K] alpha; //pseudo-counts for smoothing
    int g[N] ; // grouping factor
    int<lower=2> G; // number of groups
    int y[N]; // outcomes
}
parameters {
    simplex[K] theta[G];
}
model {
    for (j in 1:G)
        theta[j] ~ dirichlet(alpha);
    for (i in 1:N)
        y[i] ~ categorical(theta[g[i]]);
}
generated quantities {
    real score[G];
    for (k in 1:G)
        score[k] <- sum(tail(theta[k],2)) - sum(head(theta[k],7));
}
