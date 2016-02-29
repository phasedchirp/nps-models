data {
    int<lower=1> N; // number of observations
    int<lower=2> K; // number of bins
    vector[K] alpha; //pseudo-counts for smoothing
    int<lower=1,upper=K> y[N];
}
parameters {
    simplex[K] c_raw;
}
transformed parameters {
    ordered[K-1] c;
    c <- cumulative_sum(c_raw);
}
model {
    c_raw ~ dirichlet(alpha);
    for (n in 1:N)
        y[n] ~ ordered_logistic(0, c);
}
generated quantities {
    real score;
    score <- (1 - inv_logit(c[9])) - inv_logit(c[6]);
}