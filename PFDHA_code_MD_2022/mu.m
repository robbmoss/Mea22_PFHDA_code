function muout = mu(mag)
a_MD = -2.5;
b_MD = 0.415;
muout = a_MD + b_MD*mag;
end

%this is the a and b linear regression terms for P(MD|M,r)
%for all data in Mea21 the equation is  log10(MD)=0.40126*Mw-2.59939 
%for FDHI data only in Mea21 it becomes log10(MD)=0.48150*Mw-2.9305

%for complete rupture in Mea22 it becomes log10(MD)=0.415*M2-2.50
