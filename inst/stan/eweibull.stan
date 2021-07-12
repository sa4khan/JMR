functions {
matrix row_expan0(matrix X, int m){
matrix[rows(X)*m,cols(X)] Xout;
int next_row = 1;
for(i in 1:rows(X)){
Xout[next_row:(next_row + m - 1),] = rep_matrix(X[i,], m);
next_row = next_row + m;
}
return Xout;
}
matrix row_expan(matrix X, int[] m){
matrix[sum(m),cols(X)] Xout;
int next_row = 1;
for(i in 1:rows(X)){
Xout[next_row:(next_row + m[i] - 1),] = rep_matrix(X[i,], m[i]);
next_row = next_row + m[i];
}
return Xout;
}
vector row_sums(matrix X) {
vector[rows(X)] s;
for (i in 1:rows(X)) 
s[i] = sum(X[i]);
return s;
}
}
data {
int<lower=1> m;
int<lower=m> n;
int<lower=2> p1;
int<lower=1> p2;
int<lower=2> pp1;
int<lower=1> id_freq[m];
int<lower=3> qpoints;
vector<lower=0,upper=1>[m] status;
real y[n];
vector<lower=0>[m] st;
vector<lower=0,upper=1>[qpoints] weights;
matrix[m,p2] z;
matrix[m,pp1] z_rand;
matrix[n,p1] x;
matrix[n,pp1] x_rand;
matrix[m*qpoints, pp1] nodesr;
real<lower=pp1> nu;
cov_matrix[pp1] A;
vector[p1] alpha_mu;
vector<lower=0>[p1] alpha_sd;
vector[p2] beta_mu;
vector<lower=0>[p2] beta_sd;
real phi_mu;
real<lower=0> phi_sd;
real<lower=0> a0;
real<lower=0> a1;
real<lower=0> b0;
real<lower=0> b1;
real<lower=0> c0;
real<lower=0> c1;
real<lower=0> d0;
real<lower=0> d1;
}
parameters {
real<lower=0> inv_sigma2;
real<lower=0> inv_rho2;
real<lower=0> squ_kappa;
real<lower=0> squ_gam;
vector[p1] alpha;
vector[p2] beta;
real phi;
matrix[pp1,m] b_unscaled;
cov_matrix[pp1] rand_cov;
}
transformed parameters {
real<lower=0> kappa;
real<lower=0> gam;
real<lower=0> rho;
real<lower=0> sigma;
matrix[m,pp1] b;
kappa = sqrt(squ_kappa);
gam = sqrt(squ_gam);
rho = inv_sqrt(inv_rho2);
sigma = inv_sqrt(inv_sigma2);
b = (cholesky_decompose(rand_cov) * b_unscaled)';
}
model {
vector[m] pred0;
vector[m] pred1;
vector[m] psi;
vector[m] llik;
real logs0;
pred0 = z * beta;
pred1 = pred0 + phi * row_sums(z_rand .* b);
psi = (to_matrix(exp(-phi * row_sums(nodesr .* row_expan0(b,qpoints))), m, qpoints, 0) 
   * weights) .* st .* exp(-pred0) * 0.5;
for(i in 1:m){
if(status[i]){
llik[i] = log(kappa) + log(gam) + kappa * log(rho) + (kappa - 1) * log(psi[i]) 
   + (gam - 1) * log1m_exp(-pow((rho*psi[i]),kappa)) - pow((rho*psi[i]),kappa) - pred1[i];
}
else{
logs0 = log1m(exp(gam*log1m_exp(-pow((rho*psi[i]),kappa))));
llik[i] = is_inf(logs0) ? (log(gam) - pow((rho*psi[i]),kappa)) : logs0;
}
}
target += normal_lpdf(y | x * alpha + row_sums(x_rand .* row_expan(b,id_freq)), sigma);
target += llik;
target += normal_lpdf(alpha | alpha_mu, alpha_sd);
target += normal_lpdf(beta | beta_mu, beta_sd);
target += normal_lpdf(phi | phi_mu, phi_sd);
target += std_normal_lpdf(to_vector(b_unscaled));
target += inv_wishart_lpdf(rand_cov | nu, A);
target += gamma_lpdf(inv_sigma2 | a0, a1);
target += gamma_lpdf(inv_rho2 | b0, b1);
target += gamma_lpdf(squ_kappa | c0, c1);
target += gamma_lpdf(squ_gam | d0, d1);
}

