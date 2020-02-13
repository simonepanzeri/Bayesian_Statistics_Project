data
{
	int N; //observations
	vector[N] Y; //output
	real<lower=0> sigma; //variance of the residuals
	vector[3] mub;
	matrix[3,3] B;
	vector[2] G;
}

parameters {
	vector[3] beta;
}

model
{	
	real csi;
	real alpha;

	beta ~ multi_normal(mub, B);
	csi = exp(beta[1]+beta[2]*G[1]+beta[3]*G[2])/(1+exp(beta[1]+beta[2]*G[1]+beta[3]*G[2]));
	alpha = (exp(2*csi)-1)/(exp(2*csi)+1);

	// likelihood
	for(i in 2:N)
	{
		Y[i] ~ normal(alpha*Y[i-1],sigma);
	}
}
