function [index_bin]= f_CP_p_for_GLM_bins(cp, ncl)

%%%%Assign intervals of choice ratio to different choice
%%%%parameters for GLM models with multiple choice parameters. The assignment is based on a
%%%%previously calculated CP(choice ratio) profile.
%%%The same choice parameter in the GLM is assigned to intervals of choice ratio with a similar CP 

%%%Input
%%%ncl: number of choice parameters in the GLM
%%%cp: Profile CP(choice ratio)

%%%Output
%%%index_bin: index of the choice parameter assigned to each bin of the
%%%CP(choice ratio) profile

K = 50;  
cls = [];
for j = 1:K
    [ids Cs su Ds] = kmeans(cp,ncl);
    for iq = 1:ncl
        t = find(ids==iq);
        cls = [cls mean(cp(t))];
    end
end
spl = size(cp,1);
[ids Cs su Ds] = kmeans(cls',ncl);

index_bin = zeros(spl,1);
for j = 1:spl
    cos = zeros(ncl,1);
    for iq = 1:ncl
        cos(iq) = (Cs(iq)-cp(j))^2;
    end
    [te, index_bin(j)] = min(cos);
end

