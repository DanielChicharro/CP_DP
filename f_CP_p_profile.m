function [CPp, eCPp]= f_CP_p_profile(cps, ecps, psycurve, bins)

%%%%Construct a profile CP(p_cr) of choice probability as a function of
%%%%choice ratio

%%%Input 
%%% cps:       CPs for each stimulus level
%%% ecps:      standard error of CPs for each stimulus level
%%% psycurve:  psychometric function with the values of the choice ratio p_cr = p(choice D=1|stimulus level)
%%% bins:      intervals of the choice ratio p_cr for each component of the profile CP(p_cr)

%%%Output
%%%CPp: profile CP(p_cr)
%%%eCPp: standard errors of the average CP values for each bin of choice
%%%ratios

N = size(bins,1);
CPp = zeros(N,1);
eCPp = zeros(N,1);
for i = 1:N
   bin1 = bins(i,1);
   bin2 = bins(i,2);
   it1 = find(psycurve>=bin1);
   it2 = find(psycurve<bin2);
   it = intersect(it1, it2);
   icps = cps(it);
   iecps = ecps(it);
   [CPp(i), eCPp(i)]= f_weighted_average_CP(icps, iecps,0,zeros(size(icps)),zeros(size(icps)));
end