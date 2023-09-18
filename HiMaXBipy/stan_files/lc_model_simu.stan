functions {
    int is_even(int k) {
        return 2 * (k %/% 2) == k;
    }
    real median(array[] real x) {
        int K = size(x);
        if (K == 0) reject("median(x): x must not be empty");
        if (K == 1) return x[1];
        array[K] real y = sort_asc(x);
        int mid = K %/% 2;
        if (is_even(K)==1) {
            return (y[mid] + y[mid + 1]) / 2;
        } else {
            return y[mid + 1];
        }
    }
}
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
generated quantities {
    real<lower=0, upper=1> amp_dev_med;
    real<lower=0> amp_frac_med;
    real<lower=0, upper=1> amp_dev_min;
    real<lower=0> amp_frac_min;
    real<lower=0> med_rate;
    real<lower=0> max_rate;
    real<lower=0> min_rate;
    min_rate = min(sc_rate);
    med_rate = median(sc_rate);
    max_rate = max(sc_rate);
    amp_dev_med = (max_rate - med_rate) / (max_rate + med_rate);
    amp_frac_med = max_rate / med_rate;
    amp_dev_min = (max_rate - min_rate) / (max_rate + min_rate);
    amp_frac_min = max_rate / min_rate;
    array[M] real<lower=0> sc_bg_rate;
    for (m in 1:M) {
        sc_bg_rate[m] = sc_rate[m] + bg_rate[m];
    }
}