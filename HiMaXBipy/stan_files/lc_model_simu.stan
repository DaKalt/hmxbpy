data {
    int<lower=0> N;
    int<lower=0> M;
    array[M,N] int<lower=0> sc;
    array[M,N] int<lower=0> bg;
    array[M,N] real<lower=0> frac_exp;
    array[M,N] real<lower=0> dt;
    array[M,N] real<lower=0> bg_ratio;
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
    real<lower=0, upper=1> amp_dev;
    real<lower=0> amp_frac;
    real<lower=0> med_rate;
    real<lower=0> max_rate;
    med_rate = quantile(sc_rate,0.5);
    max_rate = max(sc_rate);
    amp_dev = (max_rate - med_rate) / (max_rate + med_rate);
    amp_frac = max_rate / med_rate;
    array[M] real<lower=0> sc_bg_rate;
    for (m in 1:M) {
        sc_bg_rate[m] = sc_rate[m] + bg_rate[m];
    }
}
model {
    for (m in 1:M) {
        for (n in 1:N) {
            if (dt[m,n] > 0) {
                sc[m,n] ~ poisson((sc_rate[m] + bg_rate[m]) * (dt[m,n]*frac_exp[m,n]));
                bg[m,n] ~ poisson(bg_rate[m] / bg_ratio[m,n] * (dt[m,n]*frac_exp[m,n]));
            }
        }
    }
}