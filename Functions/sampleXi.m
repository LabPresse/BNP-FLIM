function [Xi,AcceptXi] = sampleXi(Rho_IndP,Xi_Old,BNP,Chol_IndP,AcceptXi)
%sampleXi sample jointly all the means of GP priors
%
%INPUT:
%   Rho_IndP:  Learned profiles at the inducing points
%   Xi_Old:    Current values of the means of GPs
%   BNP:       Structure containing parameters used in the analysis
%   Chol_IndP: Cholesky decomposition of the correlation between inducing inducing points
%   AcceptXi:  Number of accepted proposed values so far
%
%OUTPUT:
%   Xi:        Updated GP means values
%   AcceptXi:  Updated number of accepted proposed values
%   
%Created by:
%   Mohamadreza Fazel (Presse lab 2020)
%

Sig_Xi = BNP.Sig_Xi;
Sig_Prior_Xi = BNP.Sig_Prior_Xi;
M = BNP.M;

Xi_Prop = Xi_Old+Sig_Xi*randn([1,M]);
LogLikeProp = 0;
LogLikeOld = 0;
for mm = 1:M
    LogLikeProp = LogLikeProp -0.5*sum(((log(Rho_IndP(mm).Rho)-((Xi_Prop(mm))))'/Chol_IndP).^2) ...
        - sum(log(diag(Chol_IndP))) - length(Rho_IndP(mm).Rho)*log(2*pi)/2;
    LogLikeOld = LogLikeOld -0.5*sum(((log(Rho_IndP(mm).Rho)-((Xi_Old(mm))))'/Chol_IndP).^2) ...
        - sum(log(diag(Chol_IndP))) - length(Rho_IndP(mm).Rho)*log(2*pi)/2;
end
LogLikeRatio = LogLikeProp - LogLikeOld;
LogPriorRatio = sum(log(normpdf(Xi_Prop,0,Sig_Prior_Xi))-log(normpdf(Xi_Old,0,Sig_Prior_Xi)));

if LogLikeRatio + LogPriorRatio > log(rand())
    Xi = Xi_Prop;
    AcceptXi = AcceptXi + 1;
else
    Xi = Xi_Old;
end

end
