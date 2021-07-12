data {
int<lower=1> N;
int<lower=1> ni;
int<lower=2> pp1;
int<lower=3> qpoints;
real<lower=0> st;
vector[ni*N] y;
vector[ni*N] sigma;
vector[ni] xalpha[N];
matrix[ni, pp1] z_rand;
vector[N] pred0;
vector[N] phi;
vector[N] rho;
vector[N] kappa;
row_vector<lower=0,upper=1>[qpoints] weights;
matrix[qpoints,pp1] nodes;
matrix[pp1,pp1] chol_cov_matrix[N];
}
parameters {
matrix[pp1,N] b_unscaled;
}
transformed parameters {
vector[pp1] b[N];
for(i in 1:N){
b[i] = chol_cov_matrix[i] * col(b_unscaled, i);
}
}
model {
int k = 1;
vector[ni * N] mu;
vector[N] llik;
for(i in 1:N){
mu[k:(k + ni - 1)] = xalpha[i] + z_rand * b[i];
k = k + ni;
llik[i] = lognormal_lccdf((weights * exp(-phi[i] * (nodes * b[i]))) * st 
  * exp(-pred0[i]) * 0.5 | -log(rho[i]), 1/kappa[i]);
}
target += normal_lpdf(y | mu, sigma);
target += llik;
target += std_normal_lpdf(to_vector(b_unscaled));
}

