function [mus, ps] = f_gain_model(x)

%%%Estimation of the gain model of "Partitioning neuronal variability"
%%%Goris RLT, Movshon JA, Simoncelli EP, Nature Neurosci 17, 858-865 (2014) 

%The model estimates the best fit of the spike counts to a
%negative-binomial model corresponding to Poisson responses with a gain
%fluctuation (Equation 3 Goris et al.)


%%%Input:

%%%%x: Matrix containing the spike count of the responses of a neuron. Each
%%%% different rows correspond to different stimulus levels and different
%%%% columns to different trials. Because different stimulus levels may
%%%% have a different number of trials, some value <= -1 should indicate an empty column for a certain stimulus level. 


%%%Output:

%%%mus: mean spike count for each stimulus level 

%%%ps: excitability gain factor. Assumed common for all stimulus levels.


Nc = size(x,1); %%number of stimulus levels
mus = zeros(Nc,1);%%mean spike count for each stimulus levels
s2s = zeros(Nc,1);%%%variance for each stimulus levels
nts = zeros(Nc,1);%%%number of trials for each stimulus levels
ys = [];
for j = 1:Nc
  y = x(j,x(j,:)>-1);  
  mus(j) = mean(y);
  s2s(j) = var(y);
  ys = [ys y];
  nts(j)=length(y);
end


ii = find(s2s./mus>=1);
if length(ii)>0
    x = x(ii,:); %%%selection of stimulus levels compatible with an overdispersion model
    Nc2 = size(x,1);
    mus2 = zeros(Nc2,1);
    s2s2 = zeros(Nc2,1);
    nts2 = zeros(Nc2,1);

    ys = [];
    for j = 1:Nc2
        y = x(j,x(j,:)>-1);
        mus2(j) = mean(y);
        s2s2(j) = var(y);
        ys = [ys y];
        nts2(j)=length(y);
    end

    nts22 = [0 cumsum(nts2)'];
    [ps, temp] = nbinmodelfit(ys,mus2,s2s2,nts22); %%%fit of the negative binomial model
    
else
  ps = -1;  
end
