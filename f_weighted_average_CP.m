function [CP, eCP]= f_weighted_average_CP(cps, ecps, flag, nts, p_cr)

%%%%Calculate a weighted average CP value with weights inversely
%%%%proportional to the standard error of the CP

%%%Input  
%%%cps:  array of CP values to be averaged

%%%ecps: standard error of the CP values if estimated previously

%%%flag: 0 if ecps already contains the standard errors

%%%nts: number of total trials used to calculated each CP value
%%%   needed if flag is not equal to 0

%%%p_cr: choice ratio for each CP value
%%%   needed if flag is not equal to 0


%%%Output 
%%%CP: weighted average CP
%%%eCP: standard error of the average


if abs(flag) >0
  ecps = ones(size(cps))./sqrt(12*nts.*p_cr.*(1-p_cr));  
end    

M = size(cps,1);
uw = ones(size(cps))./ecps; %% unormalized weights
nw = uw/sum(uw); %% normalized weights
CP = sum(nw.*cps);
eCP = 1/(sqrt(M)*mean(uw));
