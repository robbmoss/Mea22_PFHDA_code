function truncout = trlnrnd(mu, sigma, n)

%monte carlo sampling for a truncated lognormal distribution.  zmax is the
%maximum value where truncation occurs, here it is specified as MD = 15m.
%n random numbers between 0 and zmax are generated and thrown into the
%Inverted CDF for a truncated lognormal distribution, then pass them back.

z = unifrnd(0, 1, 1, n);

%truncout = exp(sqrt(2*sigma*sigma)*erfinv(2*z*(logncdf(log(15), mu, sigma)) - 1) + mu);

epsilon_max=5;

cdf_min=logncdf(exp(mu-epsilon_max*sigma),mu,sigma);

cdf_max=logncdf(exp(mu+epsilon_max*sigma),mu,sigma);

truncout=exp(sigma*sqrt(2)*erfinv(2*(z-0.5)*(cdf_max-cdf_min))+mu);

end

