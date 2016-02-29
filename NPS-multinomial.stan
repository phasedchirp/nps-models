data {
    int<lower=1> N; // number of observations
    int<lower=2> K; // number of bins
    vector[K] alpha; //pseudo-counts for smoothing
    int y[N]; // outcomes
}
parameters {
    simplex[K] theta;
}
model {
    theta ~ dirichlet(alpha);
    y ~ categorical(theta);
}
generated quantities {
    real score;
    score <- sum(tail(theta,2)) - sum(head(theta,7));
}