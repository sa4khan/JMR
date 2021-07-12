functions {
vector vpow(vector v, real p) {
vector[num_elements(v)] s;
for (i in 1:num_elements(v)) 
s[i] = pow(v[i],p);
return s;
}
}
data {
int<lower=1> N;
int<lower=1> ni;
int<lower=2> pp1;
int<lower=3> qpoints;
real<lower=0> st;
vector[ni*N] y;
vector[ni*N] sigma;
vector[ni] xalpha[N];
vector[qpoints] xalpha1[N];
matrix[ni, pp1] z_rand;
vector[N] pred0;
vector[N] phi;
vector[N] rho;
vector[N] kappa;
vector[N] gam;
row_vector<lower=0,upper=1>[qpoints] weights;
vector[qpoints] nodes0;
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
llik[i] = -(weights * exp((kappa[i] - 1) * log(nodes0) - log1p(vpow((gam[i] * nodes0), kappa[i]))
  + phi[i] * (xalpha1[i] + nodes * b[i]))) * st * kappa[i] * pow(rho[i],kappa[i]) * exp(pred0[i]) * 0.5;
}
target += normal_lpdf(y | mu, sigma);
target += llik;
target += std_normal_lpdf(to_vector(b_unscaled));
}

