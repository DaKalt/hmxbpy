data {
    int<lower=0> N;
    array[N] int<lower=0> sc;
    array[N] int<lower=0> bg;
    array[N] real<lower=0> frac_exp;
    array[N] real<lower=0> dt;
    array[N] real<lower=0> bg_rate;
}
parameters {
    real<lower=-3,upper=4> log_cr_sc;
    real<lower=-3,upper=4> log_cr_bg;
}
transformed parameters {
    real<lower=0> sc_rate;
    sc_rate = 10 ^ log_cr_sc;
    real<lower=0> bg_rate;
    bg_rate = 10 ^ log_cr_bg;
}
model {
    for (n in 1:N) {
        sc[n] ~ poisson((sc_rate + bg_rate) * (dt[n]*frac_exp[n]));
        bg[n] ~ poisson(bg_rate / bg_rate[n] * (dt[n]*frac_exp[n]));
    }
}