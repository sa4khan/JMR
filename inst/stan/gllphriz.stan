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
vector vpow(vector v, real p) {
vector[num_elements(v)] s;
for (i in 1:num_elements(v)) 
s[i] = pow(v[i],p);
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
matrix[m,p1] z_fixed;
matrix[n,p1] x;
matrix[n,pp1] x_rand;
vector[m*qpoints] nodes1;
matrix[m*qpoints, pp1] nodesr;
matrix[m*qpoints, p1] nodesf;
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
real<lower=0> inv_gam2;
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
gam = inv_sqrt(inv_gam2);
rho = inv_sqrt(inv_rho2);
sigma = inv_sqrt(inv_sigma2);
b = (cholesky_decompose(rand_cov) * b_unscaled)';
}
model {
vector[m] pred0;
vector[m] pred1;
pred0 = z * beta;
pred1 = pred0 + phi * (z_fixed * alpha + row_sums(z_rand .* b));
target += normal_lpdf(y | x * alpha + row_sums(x_rand .* row_expan(b,id_freq)), sigma);
target += status .* (log(kappa) + kappa * log(rho) + (kappa-1) * log(st) - log1p(vpow((gam*st),kappa)) + pred1)
   - (to_matrix(exp((kappa-1) * log(nodes1) - log1p(vpow((gam*nodes1),kappa))
   + phi * (nodesf * alpha + row_sums(nodesr .* row_expan0(b,qpoints)))), m, qpoints, 0) * weights) 
   .* st .* exp(pred0) * kappa * pow(rho,kappa) * 0.5;
target += normal_lpdf(alpha | alpha_mu, alpha_sd);
target += normal_lpdf(beta | beta_mu, beta_sd);
target += normal_lpdf(phi | phi_mu, phi_sd);
target += std_normal_lpdf(to_vector(b_unscaled));
target += inv_wishart_lpdf(rand_cov | nu, A);
target += gamma_lpdf(inv_sigma2 | a0, a1);
target += gamma_lpdf(inv_rho2 | b0, b1);
target += gamma_lpdf(squ_kappa | c0, c1);
target += gamma_lpdf(inv_gam2 | d0, d1);
}

