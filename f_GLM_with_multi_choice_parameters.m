function [bs, tdevs]= f_GLM_with_multi_choice_parameters(KK,m,grs,nu,stim,spi2,ids,choices)

%%%%Estimate a GLM model with multiple choice parameters to account for
%%%%stimulus dependent activity-choice covariations 

%%%Input
%%%KK: number of cross-validation repetitions
%%%m: polynomial order to model responses in terms of stimulus levels 
%%%grs: index of choice ratio bin for each stimulus level
%%%nu: number of choice parameters
%%%stim: stimulus level of each trial
%%%spi2: spike count of each trial
%%%ids: index of choice parameter for each choice ratio bin
%%%choices: choice for each trial, with values -1,1

%%%Output
%%%bs: estimated model parameters
%%%tdevs: likelihood of the model



t = unique(stim);
lt = length(t);

%%%create matrix of covariates and vector of the dependent variable
y = [];
X = [];
Nds = zeros(nu,2);
for k = [1:lt]
    gp = grs(k);
    gp = ids(gp);
    tt = find(stim == t(k));
    ts = spi2(tt);
    tt6 = find(isnan(ts)==0);
    tt2 = find(choices(tt(tt6))==-1);
    tt3 = find(choices(tt(tt6))==1);
    
    l2 = length(tt2);
    l3 = length(tt3);
    N = (l2+l3);
     
    Nds(gp,1) = Nds(gp,1)+l3;
    Nds(gp,2) = Nds(gp,2)+l2;
    s = spi2(tt(tt6));
    y = [y s(tt3) s(tt2)];
    
    vards = zeros(N,nu+1+m);
    vards(1:l3,gp) = 1;
    vards(l3+1:N,gp) = -1;
    for hj = 1:m
        vards(:,end-hj) = t(k)^hj;
    end
    vards(:,end) = 1;
    X = [X vards'];
end


%%%fit and cross-validate the model
tdevs = nan(KK,nu);
bs = nan(KK,nu+1+m); 

for yu = 1:KK
    yt = [];
    yf = [];
    Xt = [];
    Xf = [];
    
    nbmi = 1;
    %%%construct fitting and testing sets 
    for k = 1:nu
        N22 = min(Nds(k,:));
        N2 = round(0.8*N22); %%%test set 20%
        N3 = N22-N2;
        if N22-N2>=nbmi
            tt = find(X(k,:)==1);
            tt2 = find(X(k,:)==-1);
            re = randperm(Nds(k,1));
            re2 = randperm(Nds(k,2));
            
            yt = [yt y(tt(re(1:N3))) y(tt2(re2(1:N3)))];
            Xt = [Xt X(:,tt(re(1:N3))) X(:,tt2(re2(1:N3)))];
            yf = [yf y(tt(re(N3+1:end))) y(tt2(re2(N3+1:end)))];
            Xf = [Xf X(:,tt(re(N3+1:end))) X(:,tt2(re2(N3+1:end)))];
            
        end
    end
    
    %%%estimate the parameters and quantify the likelihood
    if isempty(yt)==0
        
        tr = find(sum(abs(Xf),2)>0);
        b = glmfit(Xf(tr,:)',yf','poisson','constant','off');
       
        ti = exp(Xt(tr,:)'*b);
        devt = poisspdf(yt,ti');
         
        ddevs = zeros(nu,1);
        for k = 1:nu
            tt = find(abs(Xt(k,:))>0);
            ddevs(k) = mean(devt(tt)); 
        end
        
        tdevs(yu,1:nu) = ddevs;
        bs(yu,tr) = b;
    end
end





