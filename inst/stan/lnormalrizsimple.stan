functions {
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
int<lower=1> which_alphafixed[p1-1];
int<lower=1> which_alphatime;
vector<lower=0,upper=1>[m] status;
real y[n];
vector<lower=0>[m] st;
matrix[m,p2] z;
matrix[m,pp1] z_rand;
matrix[m,p1] z_fixed;
matrix[n,p1] x;
matrix[m,p1-1] x_fixed;
matrix[n,pp1] x_rand;
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
}
parameters {
real<lower=0> inv_sigma2;
real<lower=0> inv_rho2;
real<lower=0> squ_kappa;
vector[p1] alpha;
vector[p2] beta;
real phi;
matrix[pp1,m] b_unscaled;
cov_matrix[pp1] rand_cov;
}
transformed parameters {
vector[p1-1] alphafixed;
real alphatime;
real<lower=0> kappa;
real<lower=0> rho;
real<lower=0> sigma;
matrix[m,pp1] b;
alphafixed = alpha[which_alphafixed];
alphatime = alpha[which_alphatime];
kappa = sqrt(squ_kappa);
rho = inv_sqrt(inv_rho2);
sigma = inv_sqrt(inv_sigma2);
b = (cholesky_decompose(rand_cov) * b_unscaled)';
}
model {
vector[m] pred0;
vector[m] pred1;
vector[m] psi;
vector[m] llik;
pred0 = z * beta;
pred1 = pred0 + phi * (z_fixed * alpha + row_sums(z_rand .* b));
psi = (exp(- pred0 - phi * (x_fixed * alphafixed + col(b, 1))) 
  .* (1 - exp(- phi * (alphatime + col(b, 2)) .* st))) ./ (phi * (alphatime + col(b, 2)));
for(i in 1:m){
llik[i] = status[i] ==1 ? (lognormal_lpdf(psi[i] | -log(rho), 1/kappa) - pred1[i]) 
   : lognormal_lccdf(psi[i] | -log(rho), 1/kappa);
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
}

