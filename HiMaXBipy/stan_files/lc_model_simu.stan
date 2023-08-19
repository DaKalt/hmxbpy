data {
    int<lower=0> N;
    int<lower=0> M;
    array[M,N] int<lower=0> sc;
    array[M,N] int<lower=0> bg;
    array[M,N] real<lower=0> frac_exp;
    array[M,N] real<lower=0> dt;
    array[M,N] real<lower=0> bg_rate;
}
parameters {
    array[M] real<lower=-3,upper=4> log_cr_sc;
    array[M] real<lower=-3,upper=4> log_cr_bg;
}
transformed parameters {
    array[M] real<lower=0> sc_rate;
    sc_rate = 10 ^ log_cr_sc;
    array[M] real<lower=0> bg_rate;
    bg_rate = 10 ^ log_cr_bg;
}
model {
    for (m in 1:M) {
        for (n in 1:N) {
            if (dt[m,n] > 0) {
                sc[m,n] ~ poisson((sc_rate[m] + bg_rate[m]) * (dt[m,n]*frac_exp[m,n]));
                bg[m,n] ~ poisson(bg_rate[m] / bg_rate[m,n] * (dt[m,n]*frac_exp[m,n]));
            }
        }
    }
}