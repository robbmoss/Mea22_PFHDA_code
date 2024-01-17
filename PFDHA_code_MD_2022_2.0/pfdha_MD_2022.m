% Probabilistic Fault Displacement Hazard Analysis (PFDHA)
% Moss et al procedure for Maximum Displacement
% Last updated 12/18/23

clear all; close all; clc;

tic
disp(' ');
disp('computing...........please wait');

%b-value for the regional seismotectonics
b_value = 0.8;

%shear modulus in dyne/cm^2
shear_modulus = 3.75*10^11;

%magnitude range for fault
min_mag = 5.0; 
max_mag = 7.5;

%Length/Width in km of fault
length = 100;
width =  15;

area = length*(1000*100)*width*(1000*100); %in cm^2

%shear wave velocity of the near surface material
vs30=700;

%location of interest along fault
%this has been descritized into 50 increments so a normalized bin range of 
%0 to 0.1 equates to xL_min=1 and xL_max=10, for a bin range of 0.4 to 0.5 
%equates to xL_min=41 and xL_max=50...
xL_min=41;
xL_max=50;

%slip rate in cm/year back calculated from IAEA rate (eq_yr) that was given
%in research request and needs to be entered below
%eq_yr=1*10^-3;

%Mo_eq=10^((3/2)*(max_mag+10.7)); %back calculating the seismic moment
%slip_rate = (Mo_eq * eq_yr) / (shear_modulus * area); %cm/yr
slip_rate= 0.5; %cm/yr

%The Truncated Exponential model is used to account for the variability in
%earthquake magnitudes 
mag = 1 : 1 : 251;

%magnitude range for this particular problem (0.01 bins)
beta = log(10)*b_value;

%number of simulations
sim=10000;

%probability density function for truncated exponential
f_m = beta * exp(-beta * (5+(2/250) * (mag - 1) - min_mag)) / (1-exp(-beta * (max_mag-min_mag)));

dm = 0.008;
dr = 0.01;

denom = 0;

for s = 1 : 1 : 250
    
denom = ( f_m(s) * 10^(1.5 * (5+(2/250) * (s - 1)) + 16.05) + f_m(s + 1) * 10^(1.5 * (5+(2/250) * (s)) + 16.05) )* dm/2 + denom;

end

N_m_min = shear_modulus*area*slip_rate / denom;


%The probability of displacement can be expressed as a function of two
%probabilities   P(D>d|M,r) = P(Slip|M,r)*P(D>d|M,r,Slip)

%The first term is the probability that fault displacement will occur given
%that an earthquake has occurred P(Slip|M,r) is modeled using the following
%function per Youngs et al. equation 4:
%P(Slip|M,r)= exp(f(x)) / (1+ exp(f(x)) where f(x)=a+bM from regression
%Youngs et al found that for all events world wide they analyzed the coefficients
%were a = -12.51 and b = 2.053
%prn = exp(a + b * (5+(2/250) * (mag-1)));
%prd = 1 + exp(a + b * (5+(2/250) * (mag-1)));
%pr_slip = prn ./ prd ;

%The logistic regression results specific to reverse events from Moss & Ross 
%2011 BSSA are in the form of the logistic function: 
%P(Slip|M,r)= 1/(1+exp(f(x)) where f(x) is the same as above.

%a=7.3;
%b=-1.03;
%pr_slip = 1./ (1 + exp(a + b * (5+(2/250) * (mag-1))));

%The logistic regression results were improved upon in Moss et al 2013 and
%an additional variable was added, VS30 of the projection of the
%rupture plane from depth.  The logistic equation is the same but the
%fuction now includes VS30 as a predictor: P(Slip|M,r,VS30)=1/(1+exp(-z))
%where z=-13.9745+2.1395*M for VS30>600m/s "stiff" materials
%      z=-6.2548+0.8308*M  for VS30<600m/s "soft" materials

if vs30>600 
pr_slip = 1./(1 + exp(-(-13.9745 + 2.1395 * (5+(2/250) * (mag-1)))));
else
pr_slip = 1./(1 + exp(-(-6.2548 + 0.8308 * (5+(2/250) * (mag-1)))));
end


%The second term in the probability of displacement defines the conditional
%distribution of the amount of fault displacement given that slip has
%occurred.  This term is analogous to an attenuation relationship and is
%constructed using empirical data.  We want a distribution that captures
%the variability of fault displacement at the site with respect to the
%entire rupture.

%We are choosing to proceed here via Monte Carlo integration, the output of
%which will be inserted into the double integral over m and r to obtain the
%final probabilities.

%There are two terms here, the probability of D/MD or D/AD at any location
%and the probability of MD or AD.  An empirical fit cumulative gamma
%distribution maps the P(D/AD) for a given x/L ratio, which is convolved
%with the P(AD) or P(MD) to arrive at the conditional probability of D.

%We now solve for the rupture displacement where the variable is treated as 
%lognormally distributed.

%m will go from 1 to 251, therefore the magnitude itself will range from
%5-7 based on the scaling relation specified below.  r will vary from 1 to
%51 and therefore vary x/L from 0 to 0.5, based on scaling relation specified
%below.

nu = zeros(28, 2);
  

probDd = zeros(51, 251, 28);
     

for m = 1 : 1 : 251 %this gives the mag bins from 5-7

    %disp('.') 

    for r = 1 : 1 : 51 %this gives the x/L bins from 0 to 0.5 with 
                       %the range for r=1:1:51 with the start value,
                       %increment, and end value
                    

        %this part creates a monte carlo for a given value of m and x
        
        %A is the Prob(MD|M) from linear regression (compared to Wells and
        %Coppersmith regression) 0.133 * 2.302 is the standard deviation
        
        %B is the Prob(D/MD|M,r) from statistical fitting of gamma or
        %weibull distributions

        A = trlnrnd(mu((5+(2/250) * (m - 1))) , 0.133*2.302, sim);
        
        B = gamrnd(a_gam(0.5*(r - 1)/50), b_gam(0.5*(r - 1)/50), 1, sim);

        B(B>=1)=1; %truncating the gamma distribution

        %This next line takes a product of the random variables A and B
        %stated above.  It is computed by simply multiplying element n of A
        %by element n of B to form combine(n).  This provides the product
        %of P(MD)*P(D/MD)

        combine = A.*B;

        
        %here a histogram is created from the data after the product of
        %random variables has occurred, with 1000 bins as the default.
        
        [n, Dbin] = hist(combine, 1000);
        
        cdf = zeros(1000,1);
        
        %CDF is created here by summing up bin weight from the left and
        %normalizing by total number of data points.
        
        cdf(1) = n(1);


        for u = 2:1:1000

            cdf(u) = n(u) + cdf(u-1);

        end

        invcdf = 1 - cdf/cdf(1000);


        %This is now the inv_CDF for a given value of x and m, but we need
        %to pull a specific value from it at d.  This wil be done by a
        %difference algorithm to grab the right value. The sensitivity is
        %set to 0.01

        %d is just a dummy index
        
        d = 1;

        
        %Here we set the range for the values of displacement we want.
        %Each value of D will yield a new rate of events exceeding D, so
        %the output of this PFDHA will yield nu(D).
        
        v = 0.01;
     
        
        D = 0.01;
        
        
        while D <= 10
      
            nu(d, 1) = D;
               
            %1000 is the last index in the CDF array.  therefore the next
            %loop goes up until 1000.  It looks at every element of the
            %Dbin array and subtracts that value from the value of D
            %specified above.  If the difference is less than the
            %sensitivity of 0.01, it pulls that value out of the array and
            %sticks it in our new discrete inverted CDF for D > D.
            
            for i = 1 : 1 : 1000

                AA = D - Dbin(i); 

                    if abs(AA) < 0.01

                        probDd(r, m, d) = invcdf(i, 1);
                        break

                    end


                % This clause in the for loop executes if the difference
                % algorithm has failed.  In other words, if no value of D
                % is found to satisfy the sensitivity algorithm it
                % Calls a spline based interpolation function--but first 
                % checks to see if D is larger than the bounds of Dbin.
                % If so it forces the probability to be zero, otherwise it
                % interpolates.


                    if i == 1000


                        if D > max(Dbin)

                            probDd(r, m, d) = 0;

                        else
                            
                            probDd(r, m, d) = interp1(Dbin, invcdf, D, 'spline');

                        end

                    end

            end
            
            
            
            %incrementing the dummy index
        
        d = d + 1;
        D = D + v;
       
        
            if D <= 0.9 && D > 0.09


                v=0.1;


            end



            if D > 0.9


                v=1;


            end
        
        
        
        end

    end



end


    %probDd now is a full Matrix.  Since we have complete arrays for
    %all three probability terms in the PFDHA integral, we now proceed to
    %calculate the rate of events per year which exceed a given value of D.

    %initializing nu to 0 before sum.

    %250 blocks here instead of 251 since there are 251 points and 250
    %actual rectangular volumes to integrate over.

    
%dd here is set to a maximum of 10 different values of D to compute nu for
%if number of D values specified above is less than 10, the remainder will
%show up as nu = 0 in the terminal.


%dd is a dummy index for a specific place in the prob(D>d) array.  it will
%go from 1 to 10 and each index corresponds to the probability for a given
%value of D from above.  So if D above went from 0.01 : 0.01 : 0.1, dd = 1
%would be D = 0.01, dd = 2 is D = 0.02, and so forth.


for dd = 1 : 1 : 28
    
    disp('.')
    
    %initializing nu as zero before sum starts
    rate = 0;
    
    for mm = 1 : 1 : 250 % magnitude range discretized into 205 increments
   
        
        for rr = xL_min : 1 : xL_max % x/L range discretized into 50 increments with
                                     % with full range of 1:1:50

            %here averaging over 4 square distanced datapoints and
            %discretely computing the integral (sum)

            rate = N_m_min*((f_m(mm) * probDd(rr, mm, dd) * pr_slip(mm) + f_m(mm) * probDd(rr+1, mm, dd) * pr_slip(mm) + f_m(mm+1) * probDd(rr+1, mm+1, dd) * pr_slip(mm+1) + f_m(mm+1) * probDd(rr, mm+1, dd) * pr_slip(mm+1))/4) * dm * dr + rate;


        end


    end

    nu(dd, 2) = rate;
    rate;

end


loglog(nu(:,1), nu(:,2))
xlabel('Displacment (m)')
ylabel('Annual Probability of Exceedance')
axis('tight');
set(gca,'FontSize',16,'FontWeight','bold')
grid on;
toc
