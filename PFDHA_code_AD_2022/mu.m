function muout = mu(mag)
a_AD = -2.87;
b_AD = 0.416;
muout = a_AD + b_AD*mag;
end

%this is the a and b linear regression terms for P(AD|M,r)
%for all data in Mea21 the equation is log10(AD)=0.39871*Mw-2.75606
%for FDHI data only the equation becomes log10(AD)=0.4395*Mw-3.0396
%Note: there is virtually no difference in these equations

%for complete rupture in Mea22          log10(AD)=0.416*Mw-2.87
