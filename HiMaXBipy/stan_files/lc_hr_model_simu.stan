data {
    int<lower=0> N;
    int<lower=0> M;
    array[M,N] int<lower=0> sc1;
    array[M,N] int<lower=0> bg1;
    array[M,N] real<lower=0> frac_exp1;
    array[M,N] int<lower=0> sc2;
    array[M,N] int<lower=0> bg2;
    array[M,N] real<lower=0> frac_exp2;
    array[M,N] real<lower=0> dt;
    array[M,N] real<lower=0> bg_ratio;
}
parameters {
    array[M] real<lower=-3,upper=4> log_cr_sc1;
    array[M] real<lower=-3,upper=4> log_cr_bg1;
    array[M] real<lower=-3,upper=4> log_cr_sc2;
    array[M] real<lower=-3,upper=4> log_cr_bg2;
}
transformed parameters {
    array[M] real<lower=0> sc_rate1;
    sc_rate1 = 10 ^ log_cr_sc1;
    array[M] real<lower=0> bg_rate1;
    bg_rate1 = 10 ^ log_cr_bg1;
    array[M] real<lower=0> sc_rate2;
    sc_rate2 = 10 ^ log_cr_sc2;
    array[M] real<lower=0> bg_rate2;
    bg_rate2 = 10 ^ log_cr_bg2;
}
model {
    for (m in 1:M) {
        for (n in 1:N) {
            if (dt[m,n] > 0) {
                sc1[m,n] ~ poisson((sc_rate1[m] + bg_rate1[m]) * (dt[m,n]*frac_exp1[m,n]));
                bg1[m,n] ~ poisson(bg_rate1[m] / bg_ratio[m,n] * (dt[m,n]*frac_exp1[m,n]));
            }
        }
    }
    for (m in 1:M) {
        for (n in 1:N) {
            if (dt[m,n] > 0) {
                sc2[m,n] ~ poisson((sc_rate2[m] + bg_rate2[m]) * (dt[m,n]*frac_exp2[m,n]));
                bg2[m,n] ~ poisson(bg_rate2[m] / bg_ratio[m,n] * (dt[m,n]*frac_exp2[m,n]));
            }
        }
    }
}
generated quantities {
    array[M] real<lower=-1,upper=1> frac;
    for (m in 1:M) {
        frac[m] = (sc_rate2[m] - sc_rate1[m]) / (sc_rate2[m] + sc_rate1[m]);
    }
    array[M] real<lower=0> sc_bg_rate1;
    for (m in 1:M) {
        sc_bg_rate1[m] = sc_rate1[m] + bg_rate1[m];
    }
    array[M] real<lower=0> sc_bg_rate2;
    for (m in 1:M) {
        sc_bg_rate2[m] = sc_rate2[m] + bg_rate2[m];
    }
}