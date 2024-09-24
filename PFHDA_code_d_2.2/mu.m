function muout = mu(mag)
a_MD = -2.5;
b_MD = 0.415;
sigma_mu = 0.148;
muout = log(10^(a_MD + b_MD*mag + sigma_mu));
end

%this is the a and b linear regression terms for P(MD|M,r)
%for all data in Mea22 the equation is  log10(MD)=0.40126*Mw-2.59939 
%for FDHI data only in Mea22 it becomes log10(MD)=0.48150*Mw-2.9305

%for complete rupture in Mea22 it becomes log10(MD)=0.415*Mw-2.50
%for uncertainty it becomes log10(MD)=0.415*Mw-2.15+(#sigma*0.148)
