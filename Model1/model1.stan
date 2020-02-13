functions {
  matrix diag_matrix (int sz)
  {
     matrix[sz,sz] M;
     for (i in 1:sz)
       for (j in 1:sz)
       {
          if(i == j) M[i,j] = 1;
          else M[i,j] = 1;
       }
     return M;
  }
  
  vector cut (int sz, vector v)
  {
     vector[sz] V;
     for (i in 1:sz)  V[i] = v[i];
     return V;
  }

  matrix cut_matrix (int sz, matrix m)
  {
     matrix[sz,sz] M;
     for (i in 1:sz)
       for (j in 1:sz)
       M[i,j] = m[i,j];
     return M;
  }

}

data {
  int<lower=1> K;
  int<lower=0> N;
  int N_int[N];
  vector[K] x1[N];
  vector[K] x2[N];
  vector[K] x3[N];
  vector[K] Y[N];
  matrix[K,K] D;
  real <lower = 0> Sigma;
  vector[3] mug;
  cov_matrix[3] Bg;
  vector[4] mub;
  cov_matrix[4] Bb;
  matrix[N,3] G;
}

parameters {
  vector[4] beta;
  vector[3] gamma;
  cov_matrix[K] Q;
}

model {
  vector[N] csi;
  vector[K] mu[N];

  beta ~ multi_normal(mub, Bb);
  gamma ~ multi_normal(mug, Bg);

  for (n in 1:N)
  { 
    mu[n] = beta[1] + beta[2] * cut(N_int[n],x1[n]) + beta[3] * cut(N_int[n],x2[n]) + beta[4] * cut(N_int[n],x3[n]);
    csi[n] = inv_logit(dot_product(gamma,G[n,]));
   }
  for (n in 1:N)
  {
    Q ~ inv_wishart(K, D*1/Sigma);
    Y[n][1:N_int[n]] ~ multi_normal_prec(mu[n], csi[n]*cut_matrix(N_int[n],Q)+(1-csi[n])*diag_matrix(N_int[n])*1/Sigma);
  }
}
