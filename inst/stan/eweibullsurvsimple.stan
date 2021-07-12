data {
int<lower=1> N;
int<lower=1> ni;
int<lower=2> pp1;
real<lower=0> st;
vector[ni*N] y;
vector[ni*N] sigma;
vector[ni] xalpha[N];
matrix[ni, pp1] z_rand;
vector[N] pred0;
vector[N] phi;
vector[N] rho;
vector[N] kappa;
vector[N] gam;
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
real psi;
real logs0;
for(i in 1:N){
mu[k:(k + ni - 1)] = xalpha[i] + z_rand * b[i];
k = k + ni;
psi = (exp(- pred0[i] - phi[i] * b[i][1]) 
 * (1 - exp(- phi[i] * b[i][2] * st))) / (phi[i] * b[i][2]);
logs0 = log1m(exp(gam[i] * log1m_exp(- pow((rho[i ] * psi), kappa[i]))));
llik[i] = is_inf(logs0) ? (log(gam[i]) - pow((rho[i] * psi), kappa[i])) : logs0;
}
target += normal_lpdf(y | mu, sigma);
target += llik;
target += std_normal_lpdf(to_vector(b_unscaled));
}

