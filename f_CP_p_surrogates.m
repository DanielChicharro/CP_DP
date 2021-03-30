function [sCPs, seCPs]= f_CP_p_surrogates(rs, choices, Nts, selstim, Nsu)


%%%Input  
%%%rs:       array with dimensions stimulus levels (rows) and maximum number of trials (columns) which contains the spike counts of all responses
%%%       choices:  array with dimension stimulus levels (rows) and maximum number of trials (columns) which contains the choice for all responses
%%%       Nts:      number of trials per stimulus level
%%%       selstim:  index of the stimulus levels selected to be used for the pool of responses to sample when constructing surrogate data sets
%%%       Nsu:      number of surrogates


%%%Output 
%%%sCPs:     CP values for all surrogates for each stimulus level
%%%       seCPs:    standard error for all surrogates for each stimulus level

z1 = []; %%pool of responses to choice D = -1
z2 = []; %%pool of responses to choice D = 1
Nsel = length(selstim);
for k = 1:Nsel
 i = selstim(k);
 irs = rs(i,1:Nts(i));
 ichoices = choices(i,1:Nts(i));
 irs1 = irs(find(ichoices==-1)); %%responses of trials with choice -1
 irs2 = irs(find(ichoices==1)); %%responses of trials with choice 1
 m1 = mean(irs1);
 m2 = mean(irs2);
 m = (m1+m2)/2;
 s1 = var(irs1);
 s2 = var(irs2);
 s = (s1+s2)/2 + ((m1-m2)^2)/4;
 izrs1 = (irs1-m)/sqrt(s); %%%normalization following Kang and Maunsell J Neurophysiol 108:3403 (2012)
 izrs2 = (irs2-m)/sqrt(s);
 z1 = [z1 izrs1];
 z2 = [z2 izrs2];
end

Nz1 = length(z1);
Nz2 = length(z2);
Nstim = size(rs,1);
sCPs = zeros(Nstim, Nsu);
seCPs = zeros(Nstim, Nsu);
for i = 1:Nstim
   ichoices = choices(i,1:Nts(i));
   n1 = length(find(ichoices==-1)); %%number of trials with choice -1
   n2 = length(find(ichoices==1));  %%number of trials with choice 1
   for j = 1:Nsu
       sz1 = zeros(n1,1);
       sz2 = zeros(n2,1);
       %%%sampling with replacement
       for k = 1:n1
          iz1 = randperm(Nz1);
          sz1(k) = z1(iz1(1));
       end
       for k = 1:n2
          iz2 = randperm(Nz2);
          sz2(k) = z2(iz2(1));
       end
       [sCPs(i,j), seCPs(i,j)]= f_CP(sz1,sz2);
   end
end

