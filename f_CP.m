function [CP, eCP]= f_CP(arr1,arr2)

%%%%Calculate the Choice Probability from the responses to two choices.


%%%Input: 
%%%arr2: contains the spike counts in the trials with choice D = 1
%%%arr1: contains the spike counts in the trials with choice D = -1

%%%Output: 
%%%CP: Choice Probability P(r_2>r_1)
%%%eCP: standard error of the CP

nb = 500; %%number of bins to estimate the probability distribution
nu1 = length(arr1);
nu2 = length(arr2);
maxr = max([max(arr1) max(arr2)]);
minr = min([min(arr1) min(arr2)]);
ib = minr-(maxr-minr)/nb:(maxr-minr)/nb:maxr+(maxr-minr)/nb;
if isempty(ib)
   CP = nan;
   eCP = nan;
else
    n1 = histc(arr1,ib)/nu1;%%probability distribution p(r|D=-1)
    n2 = histc(arr2,ib)/nu2;%%probability distribution p(r|D=1)
    tt = find(n1>0);
    nn2 = 1-cumsum(n2);
    CP = sum(n1(tt).*nn2(tt))+sum(n1(tt).*n2(tt))/2;
    eCP = 1/sqrt(12*nu1*nu2/(nu1+nu2));
end
